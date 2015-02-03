c The following can be enabled (by removing the ! before the @) to 
c speed up compilation of b2plot on IBMs where it can take a long time.
c Unfortunately the Juelich compiler complains with error message
c    (E) The STRICT/NOSTRICT option may only be specified for OPT(3).  STRICT is assumed.
c so this optimisation is disabled by default.
c 2002-09-11: logic revised
#ifdef IBMaix
#ifdef NO_OPT
@process nooptimize
#endif
#endif
      subroutine b2plot(nxd,nyd,nsd)
      use b2mod_types
      use b2mod_geo
      use b2mod_work
      use b2mod_time
      use b2mod_diag
      use b2mod_rates
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_movies
      use b2mod_sources
      use b2mod_sputter
      use b2mod_indirect
      use b2mod_elements
      use b2mod_constants
      use b2mod_residuals
      use b2mod_transport
      use b2mod_geo_corner
      use b2mod_b2plot_help
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_anomalous_transport
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
#ifdef NAGF95
      use f90_unix  	 ! IGNORE
#endif
#ifdef USE_EIRENE
      use eirmod_extrab25
      use eirmod_braeir
#endif
	implicit none
#     include <DIMENSIONS.F>
#ifndef USE_EIRENE
#     include <COUPLE/EIRDIAG.F>
#     include <COUPLE/BRAEIR.F>
#endif
      integer nxd,nyd,nsd
c
c  version : 27.07.96 20:37
c
      logical ratio
      common /t_ratio/ ratio

      real (kind=R8) ::
     *   mach(-1:nxd,-1:nyd,0:nsd-1),pr(-1:nxd,-1:nyd),
     1   fhm(-1:nxd,-1:nyd,0:1,0:nsd-1),fhp(-1:nxd,-1:nyd,0:1,0:nsd-1),
     2   fhj(-1:nxd,-1:nyd,0:1),energy_balance(-1:nxd,-1:nyd),
     3   energy_balance2(-1:nxd,-1:nyd), hy1(-1:nxd,-1:nyd),
     4   energy_balance_old(-1:nxd,-1:nyd),
     3   particle_balance(-1:nxd,-1:nyd,0:nsd-1), pb0,
     4   bas_vis_piece(-1:nxd,-1:nyd,0:nsd-1),
     *   wall_tmp(1:nwall,0:nsd-1), 
     5   trg_tmp(1:nwall,1:ndepth), 
     6   r_ref, z_ref, s_ref(-1:nyd), profile(-1:nyd),
     7   uedge_c(-1:nxd,-1:nyd), uedge_v(-1:nxd,-1:nyd,1:4),
     8   corners(-1:nxd+1,-1:nyd+1),
     9   zaf(0:nsd-1), ptf(0:nsd-1), tif, pof
      real (kind=R8) ::
     *   fetx(-1:nxd,-1:nyd), fety(-1:nxd,-1:nyd), zeff(-1:nxd,-1:nyd),
     *   fhtx(-1:nxd,-1:nyd), fhty(-1:nxd,-1:nyd)
      real (kind=R8) ::
     *   fnax(-1:nxd,-1:nyd,0:nsd-1), fnay(-1:nxd,-1:nyd,0:nsd-1),
     *   fmox(-1:nxd,-1:nyd,0:nsd-1), fmoy(-1:nxd,-1:nyd,0:nsd-1),
     *   fhmx(-1:nxd,-1:nyd,0:nsd-1), fhmy(-1:nxd,-1:nyd,0:nsd-1),
     *   fhpx(-1:nxd,-1:nyd,0:nsd-1), fhpy(-1:nxd,-1:nyd,0:nsd-1),
     *   vv(-1:nxd,-1:nyd,0:nsd-1), up(-1:nxd,-1:nyd,0:nsd-1),
     *   ww(-1:nxd,-1:nyd,0:nsd-1), ue(-1:nxd,-1:nyd),
     *   safetyfactor(-1:nyd), shear(-1:nyd), rave(-1:nyd),
     1   cvsa_cl_x(-1:nxd,-1:nyd,0:nsd-1),
     2   cvsa_cl_y(-1:nxd,-1:nyd,0:nsd-1),
     3   div_fna(-1:nxd,-1:nyd,0:nsd-1), div_fmo(-1:nxd,-1:nyd,0:nsd-1),
     3   div_fhm(-1:nxd,-1:nyd,0:nsd-1), div_fhp(-1:nxd,-1:nyd,0:nsd-1),
     3   div_fhe(-1:nxd,-1:nyd), div_fhi(-1:nxd,-1:nyd),
     3   div_fhet(-1:nxd,-1:nyd), div_fhit(-1:nxd,-1:nyd),
     4   div_fch(-1:nxd,-1:nyd), div_fhj(-1:nxd,-1:nyd),
     5   fllim0fna_x(-1:nxd,-1:nyd,0:nsd-1),
     5   fllim0fna_y(-1:nxd,-1:nyd,0:nsd-1),
     5   fllim0fhi_x(-1:nxd,-1:nyd,0:nsd-1),
     5   fllim0fhi_y(-1:nxd,-1:nyd,0:nsd-1)
#ifdef NEO
      real (kind=R8) ::                                                !cjk
     *   dn_neo(-1:nxd,-1:nyd,0:nsd), chi_neo(-1:nxd,-1:nyd,0:nsd),
     *   vn_neo(-1:nxd,-1:nyd,0:nsd), vw_neo(-1:nxd,-1:nyd,0:nsd)
      real (kind=R8) ::                                                !cjk
     *   nuab(-1:nxd,-1:nyd,0:nsd,0:nsd), nua(-1:nxd,-1:nyd,0:nsd),
     *   nustar(-1:nxd,-1:nyd,0:nsd), rhol(-1:nxd,-1:nyd,0:nsd),
     *   lnn(-1:nxd,-1:nyd,0:nsd), ltt(-1:nxd,-1:nyd,0:nsd),
     *   deltal(-1:nxd,-1:nyd,0:nsd), wb(-1:nxd,-1:nyd,0:nsd),
     *   deltaw(-1:nxd,-1:nyd,0:nsd), deltawn(-1:nxd,-1:nyd,0:nsd),
     *   deltawt(-1:nxd,-1:nyd,0:nsd), omega(-1:nxd,-1:nyd,0:nsd),
     *   deltab(-1:nxd,-1:nyd,0:nsd), deltapar(-1:nxd,-1:nyd,0:nsd),
     *   vtha(-1:nxd,-1:nyd,0:nsd)
#endif
      real (kind=R8) ::
     1   vappr_a, vappr_b, sput_frc, sput_phys, BoRiS,
     2   stick_coef, plate_kb, plate_ka, plate_thick, plate_temp,
     3   q, sput_mixed_alpha, sput_mixed_beta
      real (kind=R8) ::
     *   balance(0:nsd-1)
      real (kind=R8) ::
     *   mass,vel,vth,press
      integer ispecies(0:nsd-1), ismain, plate_model, 
     *   plate_option, boundary_namelist, neutrals_namelist, 
     *   sput_frac_flag, use_eirene

      real (kind=R8) ::
     *   rstack(1000),rpop
#ifndef NO_CDF
      real (kind=R8) ::
     *   movietime(1)
      integer :: iret, index(1), ncidcdf, varid
#endif
      integer ntstepcdf, istepcdf, lucdf
      integer istack(1000), ipop, itmp, istart, iend, itrack,
     * sput_dst, sput_src, refl_model, sput_phys_model,
     * sput_chem_model, sput_phys_col,
     * sput_dst2, sput_dst3
      character*256 sstack(1000), spop, stmp, filename, varname
      character*256 stmp1, stmp2
      integer lu, save_lu(10), ipos, save_ipos(10), lu_pos
      integer is, iss, isp, iseg, is2                            !cjk
      character*500 textline, save_textline(10)
      character*8 cval
      character*256 sval
      character*256 description
      character*22 hlp_frm
      CHARACTER FILNAM*8, H123*4, REAC*9, CRC*3
c
      real (kind=R8) ::
     *   rshift, zshift, r0, z0, a0, b0, tt0
      common /chord_shift/ rshift, zshift, r0, z0, a0, b0, tt0
c
      logical ex, cflag
      logical open_file
      integer len_of_digits
      real (kind=R8) ::
     *   qevap, qrad, 
     *   Kappa, dKappabydT, Kappa_perp, dKappa_perpbydT, Cp,
     *   mass_density, 
     *   layer_Kappa, layer_Cp, layer_emissivity
      external open_file, subini, subend, 
     *   xertst, xerrab, ipgetr, ipgeti, len_of_digits, 
     *   qevap, qrad, 
     *   Kappa, dKappabydT, Kappa_perp, dKappa_perpbydT, Cp,
     *   mass_density,
     *   layer_Kappa, layer_Cp, layer_emissivity
      external ipop, rpop, spop
      real (kind=R8) ::
     *   t1,t2
	integer, save :: background_new=0                              !srv 14.07.10
      integer it1,it2
      integer stack_req(6)
      character*256 arg
c
!!! NEED TO FIX
!!!      real (kind=R8) ::
!!!     * smfr,smpr,smth
!!!
      real (kind=R8) ::
     * b0r0,dummy,rval,sumjxa,sumjxi,suml(2),sumr(2),time 
      real (kind=R8) ::
     1 lambda, lnlam
      integer i,icell,iclv,iflag,iistack,il,imi,irstack,
     1 is1,iso,isstack,ival,ix,ix1,ix2,iy,iy1,iy2,izend,izstrt,
     2 ikind,lnblnk,nz,ixx,iyy,ixl,ixr,iyb,iyt,isel,
     3 i1,ibnd,in,ibc
      external lnblnk
#     include "ank.h"
      data sput_dst/-1/, sput_src/1/, plate_temp/300.0_R8/, ismain/1/,
     4     plate_thick/0.00_R8/, plate_model/0/, sput_frc/0.0_R8/,
     5     plate_kb/1.6e-3_R8/, plate_ka/1.3e-5_R8/, plate_option/3/,
     6     vappr_a/-40181.0_R8/, vappr_b/14.80_R8/, stick_coef/0.4_R8/,
     7     sput_phys_col/3/, sput_chem_model/0/, refl_model/1/,
     8     sput_phys_model/1/, sput_phys/0.0_R8/, lambda/-5.0_R8/,
     9     sput_dst2/-1/, sput_dst3/-1/, boundary_namelist/0/, 
     *     neutrals_namelist/0/, BoRiS/0.0_R8/, sput_frac_flag/0/,
     1     sput_mixed_alpha/1.0_R8/, sput_mixed_beta/1.0_R8/,
     2     use_eirene/0/
      save sput_dst, plate_temp, plate_thick, sput_src, sput_frc,
     1     plate_kb, plate_ka, vappr_a, vappr_b, stick_coef, ismain,
     2     plate_model, plate_option, refl_model, sput_chem_model,
     3     sput_phys_model, sput_phys_col, sput_phys, lambda, BoRiS,
     4     sput_dst2, sput_dst3, boundary_namelist, neutrals_namelist,
     5     sput_frac_flag, use_eirene

c
#ifdef USE_PXFGETENV
      integer lenval, ierror
#endif
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
*     ------------------------------------------------------------------
      ixx(ix,iy)=mod(ix+1,nx+2)-1
      iyy(ix,iy)=iy+(ix+1)/(nx+2)

c Function to check the correct parameters!
c
      call subini ('b2plot')
      call alloc_b2mod_b2plot(nxd,nyd,nsd,wklng,
     1 natm,nmol,nion,nlmpgs,nwall,ndepth)

#ifndef USE_EIRENE
      call eirene_extrab25_eirdiag_alloc(ndxp, ndyp,                     !srv 14.10.13 {
     .                     natm, nmol, nion, nlmpgs, nstra, npls)
      call eirene_extrab25_eirpbls_alloc(1, natm, nmol, nion, npls,
     &                                       nstra,nfl)
      call eirene_alloc_braeir(nxdd, nydd, nfl, 0)                       !srv 14.10.13 }
#else
      call eirene_extrab25_eirpbls_init(1, natm, nmol, nion, npls,
     .                                       nstra,nfl)
     & 
      call eirene_alloc_braeir(nxdd, nydd, nfl, 0)                       !srv 14.10.13 }
#endif

#ifdef NO_GETENV
      solpstop=' '
#else
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
      call getenv ('SOLPSTOP', solpstop)
#endif
#endif
      textan=''
      textmn=''
      textin=''
      textpl=''
      textplel=''
      textplel2=''
      textsp=''
      textwk=''
      textwwk=''
      textbl=''
      textwl=''
      texttrack=''
      do i=1,nparameters
        values(i)=0.0_R8
      end do
      time=tim
      ntstepcdf=0
      istepcdf=0
      lucdf=0
c
      ns=max(nsd,1)
      write(*,*) 'zn: ',(nint(zn(is)),is=0,ns-1)
      do is=0,ns-1
        call species(is,textpl(is))
        call species(is,textplel(is))
        call species(is,textwl(is))
      end do
      textplel(ns)=' el'
      iss=0
      il=0
      do is=0,ns-1
        if(nint(zn(is)) .ne. il) then
          il=nint(zn(is))
          if(il .le. 1) then
            imi=nint(am(is))
            if(imi .eq. 1) textsp(iss)=' '//'H'
            if(imi .eq. 2) textsp(iss)=' '//'D'
            if(imi .eq. 3) textsp(iss)=' '//'T'
          else
            textsp(iss)=' '//elements(il)
          end if
          iss=iss+1
        end if
      end do
      do is=0,ns
        do is2=0,ns
          textplel2(is,is2)=trim(textplel(is))//' -'//
     *                      trim(textplel(is2))
         enddo
      end do
      textbl(0)=' '
      textcomponents(0)=' x (poloidal)'
      textcomponents(1)=' y'
      textcomponents(2)=' z (toroidal)'
      textcomponents(3)=' (total)'
      text0123(0)=' 0'
      text0123(1)=' 1'
      text0123(2)=' 2'
      text0123(3)=' 3'
      textcorner(0)=' Bottom left'
      textcorner(1)=' Bottom right'
      textcorner(2)=' Top left'
      textcorner(3)=' Top right'
      texttrack(0)='Bulk'

!!!      if(safetyfactor(-1) .eq. 0.0_R8) call init_q

      call ipgeti ('eirene_background_new', background_new)        !srv 14.07.10
cxpb  This call is moved up from 'f31' so that we used the original cell-faced 
cxpb  quantities before their mapping to cell centers
      call ipgetr ('b2news_BoRiS', BoRiS)
      call b2xpni (nx, ny, ns, na, ni)
      nflai=0
      do is=0,ns-1
        if(is_neutral(is)) nflai=nflai+1
      enddo
      if(nflai.gt.NFL) then
        write(*,*) 'copy_background not called because nflai.gt.NFL',
     1   nflai,NFL
      elseif(nx.ge.NDXP.or.ny.ge.NDYP) then
        write(*,*)
     1   'copy_background not called because NX>=NDXP or NY>=NDYP',
     1   nx,ndxp,ny,ndyp
      else
        if(background_new.eq.0) call copy_background(nx, ny, ns, BoRiS)  !srv 14.07.10
        if(background_new.ne.0) call copy_background_(nx, ny, ns, BoRiS) !srv 14.07.10
      endif
c
      call ipgetr ('b2trcl_lambda', lambda)
      call xertst (0.0_R8.ne.lambda, 'faulty internal parameter lambda')
      call ipgeti ('b2stbr_refl_model', refl_model)
      call ipgeti ('b2stbr_plate_model', plate_model)
      if (plate_model.eq.0) then
        write(*,*) 
     .   'Using O-D time-independent plate temperature model'
      else if (plate_model.eq.1) then
        write(*,*) 'Using 1-D time-dependent plate temperature model'
      else if (plate_model.eq.2) then
        write(*,*) 'Using 2-D time-dependent plate temperature model'
      else
        call xerrab ('Invalid value for b2stbr_plate_model')
      endif
      call ipgeti ('b2stbr_plate_option', plate_option)
      call xertst (1.le.plate_option .and.
     &             3.ge.plate_option, 'faulty b2stbr_plate_option')
      call ipgeti ('b2stbr_sput_src', sput_src)
      call ipgeti ('b2stbr_sput_dst', sput_dst)
      call ipgeti ('b2stbr_sput_dst2', sput_dst2)
      call ipgeti ('b2stbr_sput_dst3', sput_dst3)
      if(sput_frc.ne.0.0_R8) then
        call xertst (0.lt.sput_src,
     1   'faulty b2stbr_sput_src')
        call xertst (sput_dst.lt.ns,
     1   'faulty b2stbr_sput_dst')
        call xertst (sput_dst2.lt.ns,
     1   'faulty b2stbr_sput_dst2')
        call xertst (sput_dst3.lt.ns,
     1   'faulty b2stbr_sput_dst3')
        call xertst (0.0_R8.le.sput_frc,
     1   'faulty b2stbr_sput_frc')
        call ipgeti ('b2stbr_sput_chem_model', sput_chem_model)
      endif
      call ipgeti ('b2stbr_sput_frac_flag', sput_frac_flag)
      call ipgetr ('b2stbr_sput_phys', sput_phys)
      call xertst (0.0_R8.le.sput_phys,
     1 'faulty input b2stbr_sput_phys')
      if(sput_phys.ne.0.0_R8) then
        call xertst (sput_dst.lt.ns,
     1   'faulty b2stbr_sput_dst')
        call ipgeti ('b2stbr_sput_phys_col', sput_phys_col)
c the TRIM sputtering files have 9 columns corresponding to
c 0 15 30 45 55 65 75 80 85 degree impact
        call xertst (0.le.sput_phys_col.and.sput_phys_col.le.9,
     1   'faulty b2stbr_sput_phys_col')
        call ipgeti ('b2stbr_sput_phys_model', sput_phys_model)
      endif
      call ipgetr ('b2stbr_sput_mixed_alpha', sput_mixed_alpha)
      call ipgetr ('b2stbr_sput_mixed_beta', sput_mixed_beta)
      call ipgetr ('b2stbr_plate_temp', plate_temp)
      call xertst (0.0_R8.lt.plate_temp, 'faulty b2stbr_plate_temp')
      call ipgetr ('b2stbr_plate_thick', plate_thick)
      call xertst (0.0_R8.le.plate_thick, 'faulty b2stbr_plate_thick')
      call ipgetr ('b2stbr_plate_kb', plate_kb)
      call xertst (0.0_R8.ne.plate_kb, 'faulty b2stbr_plate_kb')
      call ipgetr ('b2stbr_plate_ka', plate_ka)
      call xertst (0.0_R8.lt.plate_ka, 'faulty b2stbr_plate_ka')
      call ipgetr ('b2stbr_vappr_a', vappr_a)
      call ipgetr ('b2stbr_vappr_b', vappr_b)
      call ipgetr ('b2stbr_stick_coef', stick_coef)
      call xertst (0.0_R8.lt.stick_coef, 'faulty b2stbr_stick_coef')
      call ipgeti ('b2mndr_ismain', ismain)
      call xertst (0.le.ismain.and.ismain.lt.ns,
     &         'faulty b2mndr_ismain')
      call xertst (.not.is_neutral(ismain),
     &         'ismain cannot refer to a neutral species')
      call ipgeti ('b2stbc_boundary_namelist', boundary_namelist)
      call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)
      call ipgeti ('b2mndr_eirene', use_eirene)

      call initialise_potentials(ns)
      isp=0
      ispecies(0)=isp
      do is=1,ns-1
        if(zn(is) .ne. zn(is-1) .or. am(is) .ne. am(is-1)) then
          isp=isp+1
        endif
        ispecies(is)=isp
      enddo
      do is=0,ns-1
        write(*,'(i3,1x,i3,1x,f5.1,1x,f5.1,1x,f5.1,4f10.3)')
     1   is,ispecies(is),zamin(is),zamax(is),zn(is),
     2   potimin(is),potimax(is),potmin(is),potmax(is)
      end do
      call alloc_b2mod_elements(ns, zn, am)
c calculate additional fluxes (fhj, fhm, fhp)
      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ge.-1) then
            fhj(ix,iy,0)=fch(ix,iy,0)*
     1       (po(ix,iy)+po(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
          else
            fhj(ix,iy,0)=0.0_R8
          endif
          if(bottomiy(ix,iy).ge.-1) then
            fhj(ix,iy,1)=fch(ix,iy,1)*
     1       (po(ix,iy)+po(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
          else
            fhj(ix,iy,1)=0.0_R8
          endif
        enddo
      enddo
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            if(leftix(ix,iy).ge.-1) then
              fhm(ix,iy,0,is)=0.5_R8*fna(ix,iy,0,is)*
     2         (kinrgy(leftix(ix,iy),leftiy(ix,iy),is)+kinrgy(ix,iy,is))
              fhp(ix,iy,0,is)=0.5_R8*fna(ix,iy,0,is)*
     2         (rpt(leftix(ix,iy),leftiy(ix,iy),is)+rpt(ix,iy,is))*ev
            else
              fhm(ix,iy,0,is)=0.0_R8
              fhp(ix,iy,0,is)=0.0_R8
            endif
            if(bottomiy(ix,iy).ge.-1) then
              fhm(ix,iy,1,is)=0.5_R8*fna(ix,iy,1,is)*
     2         (kinrgy(bottomix(ix,iy),bottomiy(ix,iy),is)+
     1          kinrgy(ix,iy,is))
              fhp(ix,iy,1,is)=0.5_R8*fna(ix,iy,1,is)*
     2         (rpt(bottomix(ix,iy),bottomiy(ix,iy),is)+
     1          rpt(ix,iy,is))*ev
            else
              fhm(ix,iy,1,is)=0.0_R8
              fhp(ix,iy,1,is)=0.0_R8
            endif
          enddo
        enddo
      enddo
c calculate divergences of the various fluxes
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            if(leftix(ix,iy).ge.-1.and.rightix(ix,iy).le.nx.and.
     *       bottomiy(ix,iy).ge.-1.and.topiy(ix,iy).le.ny) then
              div_fna(ix,iy,is)=
     1         +fna(rightix(ix,iy),rightiy(ix,iy),0,is)-fna(ix,iy,0,is)
     2         +fna(topix(ix,iy),topiy(ix,iy),1,is)-fna(ix,iy,1,is)
              div_fmo(ix,iy,is)=
     1         +fmo(rightix(ix,iy),rightiy(ix,iy),0,is)-fmo(ix,iy,0,is)
     2         +fmo(topix(ix,iy),topiy(ix,iy),1,is)-fmo(ix,iy,1,is)
              div_fhm(ix,iy,is)=
     1         +fhm(rightix(ix,iy),rightiy(ix,iy),0,is)-fhm(ix,iy,0,is)
     2         +fhm(topix(ix,iy),topiy(ix,iy),1,is)-fhm(ix,iy,1,is)
              div_fhp(ix,iy,is)=
     1         +fhp(rightix(ix,iy),rightiy(ix,iy),0,is)-fhp(ix,iy,0,is)
     2         +fhp(topix(ix,iy),topiy(ix,iy),1,is)-fhp(ix,iy,1,is)
            else
              div_fna(ix,iy,is)=0.0_R8
              div_fmo(ix,iy,is)=0.0_R8
              div_fhm(ix,iy,is)=0.0_R8
              div_fhp(ix,iy,is)=0.0_R8
            endif
          enddo
        enddo
      enddo
      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ge.-1.and.rightix(ix,iy).le.nx.and.
     *     bottomiy(ix,iy).ge.-1.and.topiy(ix,iy).le.ny) then
            div_fhet(ix,iy)=
     1       +fhet(rightix(ix,iy),rightiy(ix,iy),0)-fhet(ix,iy,0)
     2       +fhet(topix(ix,iy),topiy(ix,iy),1)-fhet(ix,iy,1)
            div_fhit(ix,iy)=
     1       +fhit(rightix(ix,iy),rightiy(ix,iy),0)-fhit(ix,iy,0)
     2       +fhit(topix(ix,iy),topiy(ix,iy),1)-fhit(ix,iy,1)
            div_fhe(ix,iy)=
     1       +fhe(rightix(ix,iy),rightiy(ix,iy),0)-fhe(ix,iy,0)
     2       +fhe(topix(ix,iy),topiy(ix,iy),1)-fhe(ix,iy,1)
            div_fhi(ix,iy)=
     1       +fhi(rightix(ix,iy),rightiy(ix,iy),0)-fhi(ix,iy,0)
     2       +fhi(topix(ix,iy),topiy(ix,iy),1)-fhi(ix,iy,1)
            div_fch(ix,iy)=
     1       +fch(rightix(ix,iy),rightiy(ix,iy),0)-fch(ix,iy,0)
     2       +fch(topix(ix,iy),topiy(ix,iy),1)-fch(ix,iy,1)
            div_fhj(ix,iy)=0
     1       +fhj(rightix(ix,iy),rightiy(ix,iy),0)-fhj(ix,iy,0)
     2       +fhj(topix(ix,iy),topiy(ix,iy),1)-fhj(ix,iy,1)
          else
            div_fhet(ix,iy)=0.0_R8
            div_fhit(ix,iy)=0.0_R8
            div_fhe(ix,iy)=0.0_R8
            div_fhi(ix,iy)=0.0_R8
            div_fch(ix,iy)=0.0_R8
            div_fhj(ix,iy)=0.0_R8
          endif
        enddo
      enddo
c calculate the energy balance
      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ge.-1.and.rightix(ix,iy).le.nx.and.
     *     bottomiy(ix,iy).ge.-1.and.topiy(ix,iy).le.ny) then
            energy_balance(ix,iy)=-div_fhet(ix,iy)-div_fhit(ix,iy)-
     1       div_fhj(ix,iy)
            do is=0,ns-1
              energy_balance(ix,iy)=energy_balance(ix,iy)-
     1         div_fhm(ix,iy,is)*(1.0_R8-BoRiS)-div_fhp(ix,iy,is)
            enddo
c the viscous piece is still missing!
          else
            energy_balance(ix,iy)=0.0_R8
          endif
        enddo
      enddo
c additional sources
      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ge.-1.and.rightix(ix,iy).le.nx.and.
     *     bottomiy(ix,iy).ge.-1.and.topiy(ix,iy).le.ny) then
            energy_balance2(ix,iy)=
c recycling (species independent)
     1         +b2stbr_she(ix,iy)+b2stbr_shi(ix,iy)
     2         +b2stbr_sch(ix,iy)*po(ix,iy)
            do is=0,ns-1
              energy_balance2(ix,iy)=energy_balance2(ix,iy)
c radiation (electron cooling rate)
     1         -rqahe(ix,iy,is)
c recycling (species dependent)
     3         +b2stbr_smo(ix,iy,is)*ua(ix,iy,is)*(1.0_R8-BoRiS)
     4         -b2stbr_sna(ix,iy,is)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
              if(is_neutral(is)) then
                energy_balance2(ix,iy)=energy_balance2(ix,iy)
c charge exchange (neutral H/D/T)
     1           +rcxmo(ix,iy,is)*ua(ix,iy,is)*(1.0_R8-BoRiS)
     2           -rcxna(ix,iy,is)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
              else
c charge exchange from is --> is-1
                energy_balance2(ix,iy)=energy_balance2(ix,iy)
     1           -rcxmo(ix,iy,is)*ua(ix,iy,is)*(1.0_R8-BoRiS)
     2           +rcxna(ix,iy,is)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
              endif
            enddo
            do is=0,ns-2
              energy_balance2(ix,iy)=energy_balance2(ix,iy)
c ionisation from is --> is+1
     1         +rsana(ix,iy,is)*rpi(ix,iy,is+1)*ev
     2         +rsana(ix,iy,is)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
     3         -rsamo(ix,iy,is)*ua(ix,iy,is)*(1.0_R8-BoRiS)
c recombination from is+1 --> is
     4         -rrana(ix,iy,is+1)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
     5         +rramo(ix,iy,is+1)*ua(ix,iy,is)*(1.0_R8-BoRiS)
              if(.not.is_neutral(is+1)) then
                energy_balance2(ix,iy)=energy_balance2(ix,iy)
c charge exchange from is+1 --> is (receiving)
     1           -rcxna(ix,iy,is+1)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
     2           +rcxmo(ix,iy,is+1)*ua(ix,iy,is)*(1.0_R8-BoRiS)
              endif
            enddo
            do is=1,ns-1
              energy_balance2(ix,iy)=energy_balance2(ix,iy)
c recombination is --> is-1
     1         -rrana(ix,iy,is)*rpi(ix,iy,is)*ev
     2         +rrana(ix,iy,is)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
     3         -rramo(ix,iy,is)*ua(ix,iy,is)*(1.0_R8-BoRiS)
c ionisation from is-1 --> is
     4         -rsana(ix,iy,is-1)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
     5         +rsamo(ix,iy,is-1)*ua(ix,iy,is)*(1.0_R8-BoRiS)
              if(is_neutral(is-1)) then
                energy_balance2(ix,iy)=energy_balance2(ix,iy)
c charge exchange for H+/D+/T+
     1           +rcxna(ix,iy,is-1)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
     2           -rcxmo(ix,iy,is-1)*ua(ix,iy,is)*(1.0_R8-BoRiS)
              endif
            enddo
          else
            energy_balance2(ix,iy)=0.0_R8  
          endif
        enddo
      enddo
c particle balance
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            if(leftix(ix,iy).ge.-1.and.rightix(ix,iy).le.nx.and.
     *       bottomiy(ix,iy).ge.-1.and.topiy(ix,iy).le.ny) then
              pb0=-div_fna(ix,iy,is)-rsana(ix,iy,is)-rrana(ix,iy,is)+
     1         b2stbr_sna(ix,iy,is)
              if (is_neutral(is)) then
                pb0=pb0+rcxna(ix,iy,is)
              else
                pb0=pb0-rcxna(ix,iy,is)
              endif
              if(is.gt.0) then
                pb0=pb0+rsana(ix,iy,is-1)
                if(is_neutral(is-1)) then
                  pb0=pb0+rcxna(ix,iy,is-1)
                endif
              endif
              if(is.lt.ns-1) then
                pb0=pb0+rrana(ix,iy,is+1)
                if(.not.is_neutral(is+1)) then 
                  pb0=pb0+rcxna(ix,iy,is+1)
                endif
              endif
              particle_balance(ix,iy,is)=pb0
            else
              particle_balance(ix,iy,is)=0.0_R8
            endif
          enddo
        enddo
      enddo
c alternative (failed) calculation of energy balance (elegant, but flawed)
      do iy=-1,ny
        do ix=-1,nx
          energy_balance_old(ix,iy)=div_fhet(ix,iy)+div_fhit(ix,iy)-
     1     po(ix,iy)*div_fch(ix,iy)
          do is=0,ns-1
            energy_balance_old(ix,iy)=energy_balance_old(ix,iy)-
     1       ua(ix,iy,is)*div_fmo(ix,iy,is)-
     2       rpt(ix,iy,is)*ev*div_fna(ix,iy,is)+
     3       div_fhm(ix,iy,is)*(1.0_R8-BoRiS)-
     4       rqahe(ix,iy,is)
          enddo
        enddo
      enddo
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &       bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &       region(ix,iy,0).ne.0) then
              bas_vis_piece(ix,iy,is)=ua(ix,iy,is)*
     2         (cvsa_cl(rightix(ix,iy),rightiy(ix,iy),0,is)*
     4           (ua(rightix(ix,iy),rightiy(ix,iy),is)-ua(ix,iy,is))-
     5          cvsa_cl(ix,iy,0,is)*
     6           (ua(ix,iy,is)-ua(leftix(ix,iy),leftiy(ix,iy),is)))
            else
              bas_vis_piece(ix,iy,is)=0.0_R8
            endif
          enddo
        enddo
      enddo
c
      call alloc_b2mod_boundary(ns, nbcd, nnisod)
      if (boundary_namelist.eq.1)
     > call read_b2mod_boundary_namelist(nx, ny, ns, .false.)
      if(nbc.gt.nbcd) stop 'increase size in b2plot for nbcd'
      if(nniso.gt.nnisod) stop 'increase size in b2plot for nnisod'
      call alloc_b2mod_neutrals(ns, nstraid, sput_dst)
      if (neutrals_namelist.eq.1 .or. use_eirene.ne.0)
     > call read_neutrals_namelist(nx, ny, ns, .false.)
      if(nstrai.gt.nstraid) stop 'Increase DEF_NSTRA in DIMENSIONS.F'
cxpb  Bring in wall data and compute oft-needed wall quantities
      write(*,*) 'Calling alloc_b2mod_sputter'
      call alloc_b2mod_sputter(nx,ny,ns,
     1   sput_dst,sput_dst2,sput_dst3)
      call init_wall(nx,ny,ns,BoRiS,
     *   plate_model,plate_option,plate_temp,plate_thick)
      call xertst (ndepth_nml.gt.1.or.plate_model.eq.0, 
     .  'Inconsistent values of ndepth_nml and plate_model')
      call read_sputter_data(ns,
     2   sput_phys_col,sput_phys_model,refl_model)
      do i = 1, ntargsp
        texttrack(i)='Deposited '//trim(track_species(i))
      end do
      call wall_data(nx,ny,ns,BoRiS,
     *   plate_model,plate_option,plate_temp,plate_thick)
      ibnd=0
      hy1(:,:)=hy(:,:)*qz(:,:,1)
cxpb  North boundary
      do iy = ny, 0, -1
        do ix = -1, nx
          if (topiy(ix,iy).eq.ny+1 .and. leftix(ix,iy).ne.-2 .and.
     1        rightix(ix,iy).ne.nx+1 .and. region(ix,iy,0).ne.0 
     2        .and. region(ix,iy,0).ne.5) then
            ibnd=ibnd+1
            r1(ibnd)=crx(ix,iy,2)
            z1(ibnd)=cry(ix,iy,2)
            r2(ibnd)=crx(ix,iy,3)
            z2(ibnd)=cry(ix,iy,3)
            wall_heatflux(ibnd)=max(0.0_R8,fhi(ix,iy,1))/gs(ix,iy,1)
            wall_area(ibnd)=gs(ix,iy,1)
*    ..compute cell-faced quantities
            do is = 0, ns-1
               zaf(is) = 
     .             (rza(ix,iy,is)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &              rza(bottomix(ix,iy),bottomiy(ix,iy),is)*hy1(ix,iy))/
     &             (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
               ptf(is) = 
     .             (rpt(ix,iy,is)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &              rpt(bottomix(ix,iy),bottomiy(ix,iy),is)*hy1(ix,iy))/
     &             (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
            enddo
            tif = (ti(ix,iy)*hy(bottomix(ix,iy),bottomiy(ix,iy))+
     &        ti(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &       (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
            pof = (po(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &        po(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &       (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
            do is = 0, ns-1
              incident_energy(ibnd,is)=
     .         max(0.0_R8,tif/ev+pof*zaf(is)) ! xpb: Incident energy (eV)
              wall_flux(ibnd,is)=max(0.0_R8,fna(ix,iy,1,is))
              wall_heatflux(ibnd)=wall_heatflux(ibnd)+
     .         max(0.0_R8,fna(ix,iy,1,is))/gs(ix,iy,1)*
     &          (kinrgy(ix,iy,is)*(1.0_R8-BoRiS) + ptf(is)*ev)
            enddo
            wall_heatflux(ibnd)=wall_heatflux(ibnd)+
     .       max(0.0_R8,fhe(ix,iy,1))/gs(ix,iy,1) 
          endif
        enddo
      enddo
cxpb  East boundary
      do ix = nx, 0, -1
        do iy = ny, -1, -1
          if (rightix(ix,iy).eq.nx+1 .and. bottomiy(ix,iy).ne.-2
     1         .and. topiy(ix,iy).ne.ny+1 .and. region(ix,iy,0).ne.0)
     >         then
            ibnd=ibnd+1
            r1(ibnd)=crx(ix,iy,3)
            z1(ibnd)=cry(ix,iy,3)
            r2(ibnd)=crx(ix,iy,1)
            z2(ibnd)=cry(ix,iy,1)
            wall_heatflux(ibnd)=max(0.0_R8,fhi(ix,iy,0))/gs(ix,iy,0)
            wall_area(ibnd)=gs(ix,iy,0)
*    ..compute cell-faced quantities
            do is = 0, ns-1
               zaf(is) = 
     .             (rza(ix,iy,is)*hx(leftix(ix,iy),leftiy(ix,iy))+
     &              rza(leftix(ix,iy),leftiy(ix,iy),is)*hx(ix,iy))/
     &             (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
               ptf(is) = 
     .             (rpt(ix,iy,is)*hx(leftix(ix,iy),leftiy(ix,iy))+
     &              rpt(leftix(ix,iy),leftiy(ix,iy),is)*hx(ix,iy))/
     &             (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
            enddo
            tif = (ti(ix,iy)*hx(leftix(ix,iy),leftiy(ix,iy))+
     &        ti(leftix(ix,iy),leftiy(ix,iy))*hx(ix,iy))/
     &       (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
            pof = (po(ix,iy)*hx(leftix(ix,iy),leftiy(ix,iy))+
     &        po(leftix(ix,iy),leftiy(ix,iy))*hx(ix,iy))/
     &       (hy1(leftix(ix,iy),leftiy(ix,iy))+hy1(ix,iy))
            do is = 0, ns-1
              incident_energy(ibnd,is)=
     .         max(0.0_R8,tif/ev+pof*zaf(is)) ! xpb: Incident energy (eV)
              wall_flux(ibnd,is)=max(0.0_R8,fna(ix,iy,0,is))
              wall_heatflux(ibnd)=wall_heatflux(ibnd)+
     .         max(0.0_R8,fna(ix,iy,0,is))/gs(ix,iy,0)*
     &         (kinrgy(ix,iy,is)*(1.0_R8-BoRiS) + ptf(is)*ev)
            enddo
            wall_heatflux(ibnd)=wall_heatflux(ibnd)+
     .       max(0.0_R8,fhe(ix,iy,0))/gs(ix,iy,0)  
          endif
        enddo
      enddo
cxpb  South Boundary
      do iy = -1, ny-1
        do ix = nx, -1, -1
          if(bottomiy(ix,iy).eq.-2 .and. 
     1           leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and. 
     2           mod(region(ix,iy,0),4).ne.1 .and. region(ix,iy,0).ne.0) 
     >            then
            ibnd=ibnd+1
            r1(ibnd)=crx(ix,iy,1)
            z1(ibnd)=cry(ix,iy,1)
            r2(ibnd)=crx(ix,iy,0)
            z2(ibnd)=cry(ix,iy,0)
*    ..compute cell-faced quantities
            do is = 0, ns-1
               zaf(is) = 
     .             (rza(ix,iy,is)*hy1(topix(ix,iy),topiy(ix,iy))+
     &              rza(topix(ix,iy),topiy(ix,iy),is)*hy1(ix,iy))/
     &             (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
               ptf(is) = 
     .             (rpt(ix,iy,is)*hy1(topix(ix,iy),topiy(ix,iy))+
     &              rpt(topix(ix,iy),topiy(ix,iy),is)*hy1(ix,iy))/
     &             (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
            enddo
            tif = (ti(ix,iy)*hy(topix(ix,iy),topiy(ix,iy))+
     &        ti(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &       (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
            pof = (po(ix,iy)*hy1(topix(ix,iy),topiy(ix,iy))+
     &        po(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &       (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
            wall_heatflux(ibnd)=
     .       max(0.0_R8,-fhi(topix(ix,iy),topiy(ix,iy),1))/
     &             gs(topix(ix,iy),topiy(ix,iy),1)
            wall_area(ibnd)=gs(topix(ix,iy),topiy(ix,iy),1)
            do is = 0, ns-1
              incident_energy(ibnd,is)=
     .         max(0.0_R8,tif/ev+pof*zaf(is)) ! xpb: Incident energy (eV)
              wall_flux(ibnd,is)=max(0.0_R8,
     .         -fna(topix(ix,iy),topiy(ix,iy),1,is))
              wall_heatflux(ibnd)=wall_heatflux(ibnd) + max(0.0_R8,
     &                  -fna(topix(ix,iy),topiy(ix,iy),1,is))/
     &                    gs(topix(ix,iy),topiy(ix,iy),1)*
     &              (kinrgy(ix,iy,is)*(1.0_R8-BoRiS) + ptf(is)*ev)
            enddo
            wall_heatflux(ibnd)=wall_heatflux(ibnd)+
     .         max(0.0_R8,-fhe(topix(ix,iy),topiy(ix,iy),1))/
     &             gs(topix(ix,iy),topiy(ix,iy),1)
          endif
        enddo
      enddo
cxpb  Western boundary
      do ix = -1, nx-1
        do iy = -1, ny
          if (leftix(ix,iy).eq.-2 .and. bottomiy(ix,iy).ne.-2 .and.
     1         topiy(ix,iy).ne.ny+1 .and. region(ix,iy,0).ne.0) then
            ibnd=ibnd+1
            r1(ibnd)=crx(ix,iy,0)
            z1(ibnd)=cry(ix,iy,0)
            r2(ibnd)=crx(ix,iy,2)
            z2(ibnd)=cry(ix,iy,2)
            wall_heatflux(ibnd)=max(0.0_R8,
     &             -fhi(rightix(ix,iy),rightiy(ix,iy),0))/
     &             gs(rightix(ix,iy),rightiy(ix,iy),0)
            wall_area(ibnd)=gs(rightix(ix,iy),rightiy(ix,iy),0)
*    ..compute cell-faced quantities
            do is = 0, ns-1
               zaf(is) = 
     .             (rza(ix,iy,is)*hx(rightix(ix,iy),rightiy(ix,iy))+
     &              rza(rightix(ix,iy),rightiy(ix,iy),is)*hx(ix,iy))/
     &             (hx(rightix(ix,iy),rightiy(ix,iy))+hx(ix,iy))
               ptf(is) = 
     .             (rpt(ix,iy,is)*hx(rightix(ix,iy),rightiy(ix,iy))+
     &              rpt(rightix(ix,iy),rightiy(ix,iy),is)*hx(ix,iy))/
     &             (hx(rightix(ix,iy),rightiy(ix,iy))+hx(ix,iy))
            enddo
            tif = (ti(ix,iy)*hx(rightix(ix,iy),rightiy(ix,iy))+
     &        ti(rightix(ix,iy),rightiy(ix,iy))*hx(ix,iy))/
     &       (hx(rightix(ix,iy),rightiy(ix,iy))+hx(ix,iy))
            pof = (po(ix,iy)*hx(rightix(ix,iy),rightiy(ix,iy))+
     &        po(rightix(ix,iy),rightiy(ix,iy))*hx(ix,iy))/
     &       (hy1(rightix(ix,iy),rightiy(ix,iy))+hy1(ix,iy))
            do is = 0, ns-1
              incident_energy(ibnd,is)=
     .         max(0.0_R8,tif/ev+pof*zaf(is)) ! xpb: Incident energy (eV)
              wall_flux(ibnd,is)=max(0.0_R8,
     .         -fna(rightix(ix,iy),rightiy(ix,iy),0,is))
              wall_heatflux(ibnd)=wall_heatflux(ibnd)+
     &             max(0.0_R8,-fna(rightix(ix,iy),rightiy(ix,iy),0,is))/
     &             gs(rightix(ix,iy),rightiy(ix,iy),0)*
     &             (kinrgy(ix,iy,is)*(1.0_R8-BoRiS) + ptf(is)*ev)
            enddo
            wall_heatflux(ibnd)=wall_heatflux(ibnd)+
     .       max(0.0_R8,-fhe(rightix(ix,iy),rightiy(ix,iy),0))/
     &             gs(rightix(ix,iy),rightiy(ix,iy),0)  
          endif
        enddo
      enddo
      ndata = nbnd
      do ibnd = 1, nbnd
        wall_temp(ibnd) = target_temp(ibnd,1)
      enddo
      target_stack(1:nbnd,1:ndepth_nml,1) = 
     1  target_temp(1:nbnd,1:ndepth_nml)
      target_label(1)='Target temperature (K)'
      call set_targetbarlabel(target_label(1))
      if (ndepth_nml.gt.1) then
        do ibnd = 1, nbnd
          target_vol(ibnd,1) = target_width(ibnd,1)*
     1     ((target_depth(ibnd,2)-target_depth(ibnd,1))/2.0_R8)*
     2     pi*(targetcr(ibnd,1,1)+targetcr(ibnd,1,2))
          do is = 2, ndepth_nml-1
            target_vol(ibnd,is) = target_width(ibnd,is)*
     1       (target_depth(ibnd,is)-target_depth(ibnd,is-1))*
     2       pi*(targetcr(ibnd,is,1)+targetcr(ibnd,is,2))   
          enddo
          target_vol(ibnd,ndepth_nml) = target_width(ibnd,ndepth_nml)*
     1     ((target_depth(ibnd,ndepth_nml)-
     .       target_depth(ibnd,ndepth_nml-1))
     2      /2.0_R8)*
     3     pi*(targetcr(ibnd,ndepth_nml,1)+targetcr(ibnd,ndepth_nml,2))
        enddo
      else
        target_vol(1:nbnd,1) = target_width(1:nbnd,1)
      endif
c
      call mapx(nx,ny,1,gs(:,:,0),work)
      call mapy(nx,ny,1,gs(:,:,1),work)
!!!      call mapx(nx,ny,ns,ua,work)
      do is=0,ns-1
        call mapx(nx,ny,1,fna(:,:,0,is),work)
        call mapy(nx,ny,1,fna(:,:,1,is),work)
        call mapx(nx,ny,1,fna_32(:,:,0,is),work)
        call mapy(nx,ny,1,fna_32(:,:,1,is),work)
        call mapx(nx,ny,1,fna_52(:,:,0,is),work)
        call mapy(nx,ny,1,fna_52(:,:,1,is),work)
        call mapx(nx,ny,1,fna_nodrift(:,:,0,is),work)
        call mapy(nx,ny,1,fna_nodrift(:,:,1,is),work)
        call mapx(nx,ny,1,fmo(:,:,0,is),work)
        call mapy(nx,ny,1,fmo(:,:,1,is),work)
        call mapx(nx,ny,1,cvsa_cl(:,:,0,is),work)
        call mapy(nx,ny,1,cvsa_cl(:,:,1,is),work)
        call mapx(nx,ny,1,fllim0fna(:,:,0,is),work)
        call mapy(nx,ny,1,fllim0fna(:,:,1,is),work)
        call mapx(nx,ny,1,fllim0fhi(:,:,0,is),work)
        call mapy(nx,ny,1,fllim0fhi(:,:,1,is),work)
        call mapx(nx,ny,1,fhm(:,:,0,is),work)
        call mapy(nx,ny,1,fhm(:,:,1,is),work)
        call mapx(nx,ny,1,fhp(:,:,0,is),work)
        call mapy(nx,ny,1,fhp(:,:,1,is),work)
      enddo
      call mapx(nx,ny,1,fhe(:,:,0),work)
      call mapy(nx,ny,1,fhe(:,:,1),work)
      call mapx(nx,ny,1,fhi(:,:,0),work)
      call mapy(nx,ny,1,fhi(:,:,1),work)
      call mapx(nx,ny,1,fhet(:,:,0),work)
      call mapy(nx,ny,1,fhet(:,:,1),work)
      call mapx(nx,ny,1,fhit(:,:,0),work)
      call mapy(nx,ny,1,fhit(:,:,1),work)
      call mapx(nx,ny,1,fhj(:,:,0),work)
      call mapy(nx,ny,1,fhj(:,:,1),work)
      call mapx(nx,ny,1,fch(:,:,0),work)
      call mapy(nx,ny,1,fch(:,:,1),work)
      call mapx(nx,ny,1,fchdia(:,:,0),work)  	 ! cd
      call mapy(nx,ny,1,fchdia(:,:,1),work)  	 ! dd
      call mapx(nx,ny,1,fch_32(:,:,0),work)
      call mapy(nx,ny,1,fch_32(:,:,1),work)
      call mapx(nx,ny,1,fch_52(:,:,0),work)
      call mapy(nx,ny,1,fch_52(:,:,1),work)
      call mapx(nx,ny,1,fchvispar(:,:,0),work)
      call mapy(nx,ny,1,fchvispar(:,:,1),work)
      call mapx(nx,ny,1,fchvisper(:,:,0),work)
      call mapy(nx,ny,1,fchvisper(:,:,1),work)
      call mapx(nx,ny,1,fchin(:,:,0),work)
      call mapy(nx,ny,1,fchin(:,:,1),work)
      call mapx(nx,ny,1,fch_p(:,:,0),work)           !srv 10.06.09
      call mapy(nx,ny,1,fch_p(:,:,1),work)           !srv 10.06.09
!!!      call mapx(nx,ny,ns,smpr,work)
!!!      call mapx(nx,ny,ns,smth,work)
!!!      call mapx(nx,ny,ns,smfr,work)
c
      inquire(file='B2SDIAG',exist=ex)
#ifndef HPUX
      inquire(file='fort.44',exist=ex)
#else
      inquire(file='ftn44',exist=ex)
#endif
      if(ex) then
        call ngread(ismain,iflag)
        textwl(ns:ns+natmi-1)=textan(0:natmi-1)
        textwl(ns+natmi:ns+natmi+nmoli-1)=textmn(0:nmoli-1)
        textwl(ns+natmi+nmoli:ns+natmi+nmoli+nioni-1)=textin(0:nioni-1)
        nsw=ns+natmi+nmoli+nioni
      else
        natmi=1
        textan(0)='ZERO'
        nmoli=1
        textmn(0)='ZERO'
        nioni=1
        textin(0)='ZERO'
        nsw=ns
        call eirene_mc_init
        call alloc_b2mod_diag(nx,ny,ns,nnreg(0),nstraid,natm)
      endif
c
      nh=0.0_R8
      do is=0,ns-1
        stmp1=trim(textpl(is))
        call lower_case(stmp1)                                !srv 14.10.13
        if(is_neutral(is) .and. nint(zn(is)).eq.1) then
          nh(-1:nx,-1:ny)=nh(-1:nx,-1:ny)+na(-1:nx,-1:ny,is)
          do iss=1,natmi
            stmp2=trim(textan(iss-1))//'0'
            call lower_case(stmp2)                                !srv 14.10.13
            if(stmp1.eq.stmp2) then
              nh(-1:nx,-1:ny)=nh(-1:nx,-1:ny)+dab2(0:nx+1,0:ny+1,iss,1)
            endif
          enddo
        endif
      enddo
c
c this just sets initialised_plotting to false --- the actual plotting
c initialisation is done on the first call to init_page so that color maps
c and page sizes can be set before plotting is started
c
      call initialise_plotting_flags(time)
c
      nz=ns
c
      do ix=-1,nx
        XX(IX,2) = ix
        XX(IX,1) = IX-0.5_R8
      end do
      do iy=-1,ny
        YY(IY,2) = IY
        YY(IY,1) = IY-0.5_R8
      end do
      do ibnd=1,nbnd
        targetx(ibnd,2) = real(ibnd)
        targetx(ibnd,1) = real(ibnd)-0.5_R8
      end do
      do i=1,ndepth_nml
        targety(i,2) = real(i)
        targety(i,1) = real(i)-0.5_R8
      end do
      
      pr=ne*te+sum(na,3)*ti
      cr=0.25_R8*sum(crx,3)
      cz=0.25_R8*sum(cry,3)
      targetr=0.5_R8*sum(targetcr,3)
      targetz=0.5_R8*sum(targetcz,3)
      icell=0
      r0=0.0_R8
      z0=0.0_R8
      do ix = -1, nx
        if (mod(region(ix,-1,0),4).eq.1) then
          icell=icell+1
          r0=r0+cr(ix,-1)
          z0=z0+cz(ix,-1)
        endif
      enddo
      if (icell.gt.0) then
        r0=r0/icell
        z0=z0/icell
      else
        r0=cr(-1,-1)
        z0=cz(-1,-1)
      endif
      values(1)=minval(xx)                                  !xmin
      values(2)=maxval(xx)+0.5_R8                           !xmax
      values(3)=minval(yy)                                  !ymin
      values(4)=maxval(yy)+0.5_R8                           !ymax
      values(5)=minval(crx)                                 !rmin
      values(6)=maxval(crx)                                 !rmax
      values(7)=minval(cry)                                 !zmin
      values(8)=maxval(cry)                                 !zmax
      values(11)=init_value                           	 !omin
      values(12)=init_value                            	 !omax
      values(13)=0.5_R8                                     !xref
      values(16)=5.0e-4_R8                                  !traj: plottime
      values(17)=0.0_R8                                     !rshift
      values(18)=0.0_R8                                     !zshift
      values(19)=0.0_R8                                     !r0
      values(20)=0.0_R8                                     !z0
      values(21)=0.0_R8                                     !a0
      values(22)=0.0_R8                                     !b0
      values(23)=0.0_R8                                     !tt0
      values(24)=minval(targetx(1:nbnd,:))                  !xmin target
      values(25)=maxval(targetx(1:nbnd,:))+0.5_R8           !xmax target
      values(26)=minval(targety(1:ndepth_nml,:))            !ymin target
      values(27)=maxval(targety(1:ndepth_nml,:))+0.5_R8     !ymax target
      values(28)=minval(targetcr(1:nbnd,:,:))               !rmin target
      values(29)=maxval(targetcr(1:nbnd,:,:))               !rmax target
      values(30)=minval(targetcz(1:nbnd,:,:))               !zmin target
      values(31)=maxval(targetcz(1:nbnd,:,:))               !zmax target
      values(32)=minval(target_temp(1:nbnd,:))              !fmin target
      values(33)=maxval(target_temp(1:nbnd,:))              !fmax target
      values(34)=init_value                          	 !omin target
      values(35)=init_value                           	 !omax target
      values(36)=init_value                           	 !smin wall
      values(37)=init_value                            	 !smax wall
      expfile=' '
      wallfile=' '
c
      call b2page(1,1)
      color_set=.true.
      call contourplot(.true.)     ! draw contour plot by default
      call nofill(.false.)	 ! fill contours by default
      call ipolate(.false.)        ! do not use splines untless wanted
      call clvl(.true.)   	 ! contour level lines plotted
      call labelon(.true.)	 ! contour labels drawn
      call gridon(.false.)	 ! overlayed mesh off by default
      call vesselon(.false.)       ! same for vessel
      call separatrixon(.false.)   ! same for separatrix
      call set_midplane_mapping(.false.) ! do not map to midplane by default
      call wallon(.false.,0,values(14),values(15),1)
      call targeton(.false.,values(32),values(33),1)
      call set_partition_colors(0.0_R8)
      plotting_wall=.false.
      b2plot_debug=0
      sepwidth=1.0_R8     	 ! standard separatrix line width
      vwidth=1.0_R8       	 ! standard vessel line width
      owidth=1.0_R8       	 ! standard outline line width
      ewidth=1.0_R8       	 ! standard experimental data line width
      arroww=1.0_R8       	 ! standard arrow line width
      lwidth=1.0_R8       	 ! linewidth
      wwidth=10.0_R8      	 ! wall linewidth
      zwidth=1.0_R8              ! axis linewidth
      lfsz=0                     ! line label font size
      typeofplot=3        	 ! physical and computational
      xparallel=.false.   	 ! default to f.y using poloidal distance
      call set_fill_option(0)
      linlog=1            	 ! linear plot
      wlinlog=1           	 ! linear plot
      tlinlog=0           	 ! log plot for trajectories
      trglinlog=1         	 ! linear plot for targets
      call set_square     	 ! square plot
      writing_to_file=.false.
      call initialise_plotting_0
      moviename=' '
      movie_skip=1
      average=1
      use_compile=.true.
      acol=1
      bcol=1
      dcol=1
      ecol=1
      gcol=1
      hcol=1
      mcol=1
      ocol=1
      scol=1
      tcol=1
      vcol=1
      wcol=1
      xcol=1
      zcol=1
      set_resignore_regions_to_zero=.false.
      ue(:,:)=-fch_p(:,:,0)/gs(:,:,0)/qc(:,:)/qe       !srv 10.06.09
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            vv(ix,iy,is)=fna(ix,iy,1,is)/na(ix,iy,is)/gs(ix,iy,1)
            up(ix,iy,is)=
     .       fna(ix,iy,0,is)/na(ix,iy,is)/gs(ix,iy,0)/qc(ix,iy)
            ue(ix,iy)=ue(ix,iy)+rza(ix,iy,is)*ua(ix,iy,is)*na(ix,iy,is)
            fnax(ix,iy,is)=fna(ix,iy,0,is)
            fmox(ix,iy,is)=fmo(ix,iy,0,is)
            fhmx(ix,iy,is)=fhm(ix,iy,0,is)
            fhpx(ix,iy,is)=fhp(ix,iy,0,is)
            arrowx(ix,iy,is)=fna(ix,iy,0,is)
            fnay(ix,iy,is)=fna(ix,iy,1,is)
            fmoy(ix,iy,is)=fmo(ix,iy,1,is)
            fhmy(ix,iy,is)=fhm(ix,iy,1,is)
            fhpy(ix,iy,is)=fhp(ix,iy,1,is)
            arrowy(ix,iy,is)=fna(ix,iy,1,is)
            uutraj(ix,iy,is)=ua(ix,iy,is)
            vvtraj(ix,iy,is)=vv(ix,iy,is)
            cvsa_cl_x(ix,iy,is)=cvsa_cl(ix,iy,0,is)
            cvsa_cl_y(ix,iy,is)=cvsa_cl(ix,iy,1,is)
            fllim0fna_x(ix,iy,is)=fllim0fna(ix,iy,0,is)
            fllim0fna_y(ix,iy,is)=fllim0fna(ix,iy,1,is)
            fllim0fhi_x(ix,iy,is)=fllim0fhi(ix,iy,0,is)
            fllim0fhi_y(ix,iy,is)=fllim0fhi(ix,iy,1,is)
          enddo
        enddo
      enddo
      ue(:,:)=ue(:,:)/ne(:,:)
      do iy=-1,ny
        rave(iy)=0.0_R8
        t1=0.0_R8
        t2=0.0_R8
        icell=0
        q=0.0_R8
        do ix=-1,nx
          if(mod(region(ix,iy,0),4).eq.1) then
            q=q+abs(bb(ix,iy,3)/bb(ix,iy,0))/cr(ix,iy)*hx(ix,iy)
            t1=t1+hx(ix,iy)*bb(ix,iy,3)/abs(bb(ix,iy,0))
            t2=t2+hx(ix,iy)
            rave(iy)=rave(iy)+sqrt((cr(ix,iy)-r0)**2+
     &                             (cz(ix,iy)-z0)**2)
            icell=icell+1
          endif
        enddo
        q=q/2.0_R8/pi
        if (icell.gt.0) then
          safetyfactor(iy)=t1/t2
          write(*,*) 'safetyfactor,q: ',safetyfactor(iy),q
          safetyfactor(iy)=q    ! I think this is better [DPC, 2004-09-24]
          rave(iy)=rave(iy)/icell
        else
          safetyfactor(iy)=0.0_R8
        endif
      enddo
      do iy=-1,ny
        if (iy.eq.-1) then
          shear(iy)=(safetyfactor(iy+1)-safetyfactor(iy))/
     &     (rave(iy+1)-rave(iy))
        elseif (iy.lt.jsep) then
          shear(iy)=(safetyfactor(iy+1)-safetyfactor(iy-1))/
     &     (rave(iy+1)-rave(iy-1))
        elseif (iy.eq.jsep) then
          shear(iy)=(safetyfactor(iy)-safetyfactor(iy-1))/
     &     (rave(iy)-rave(iy-1))
        else
          shear(iy)=0.0_R8
        endif
      enddo

c calculate neoclassical diagnostics
#ifdef NEO      
      call neo_diagnostics(nx, ny, ns, dn_neo, chi_neo, 
     &   vn_neo, vw_neo)
      call neo_validity(nx, ny, ns, jsep, jxa, cr, safetyfactor,
     &   nuab, nua, nustar, rhol, lnn, ltt, deltal, wb, deltaw, 
     &   deltawn, deltawt, omega, deltab, deltapar, vtha)
#endif
c
      lu=5
      lu_pos=0
      call readline(textline,lu,ipos)
c
      irstack=0
      iistack=0
      isstack=0
      iwrk=0
      iwwrk=0
      itwrk=1

   10 continue
      if(b2plot_debug.ge.5)
     1 write(*,'(a,6i4)') 'STACK i,r,s,m,w,t = ',
     . iistack,irstack,isstack,iwrk,iwwrk,itwrk
      if(b2plot_debug.ge.12)
     1 write(*,*) 'ISTACK: ',iistack,(istack(i),i=1,iistack)
      if(b2plot_debug.ge.12)
     1 write(*,*) 'RSTACK: ',irstack,(rstack(i),i=1,irstack)
      call symbol(textline,lu,ipos,ikind,ival,rval,sval,cval)
      if(ikind .eq. 1) then
        call ipush(istack,iistack,ival)
        goto 10
      else if(ikind .eq. 2) then
        call rpush(rstack,irstack,rval)
        goto 10
      else if(ikind .eq. 3) then
        call spush(sstack,isstack,sval)
        if(b2plot_debug.ge.6)
     1   write(*,*) 'Read string: ',sval(1:lnblnk(sval))
        goto 10
      end if
      call lower_case(cval)                                !srv 14.10.13

      cflag=.false.
      if (.not.cflag) call select_matrix
      if (.not.cflag) call select_stack_manipulation
      if (.not.cflag) call select_parameters
      if (.not.cflag) call select_command
      if (.not.cflag) call select_wall_stack
      if (.not.cflag) call select_wall_command
      if (.not.cflag) call select_target_stack
      if (.not.cflag) call select_constants
      if (.not.cflag) call select_flvalue

      if (.not.cflag) then
        if (cval.eq.'eof') then    ! end of file
          if (lu.eq.5) then
            call terminate_plotting
            write(*,'(a,6i4)') 'STACK i,r,s,m,w,t = ',
     .       iistack,irstack,isstack,iwrk,iwwrk,itwrk
            call subend ()
            return
          end if
#ifndef NO_CDF
          if(ntstepcdf.gt.0.and.istepcdf.lt.ntstepcdf
     1     .and.lucdf.eq.lu) then
            rewind(lu)
            istepcdf=istepcdf+1
            call get_from_movie_cdf(ncidcdf,istepcdf)
            rewind(lu)
          elseif(ntstepcdf.gt.0.and.istepcdf.eq.ntstepcdf
     1     .and.lucdf.eq.lu) then
            iret = nf_close (ncidcdf)
            call check_cdf_status(iret)
            write(*,*) 'Finished with cmdcdf'
            close(lu)
            lu=save_lu(lu_pos)
            ipos=save_ipos(lu_pos)
            textline=save_textline(lu_pos)
            lu_pos=lu_pos-1
            write(*,*) 'returning to lu ',lu
          else
            close(lu)
            lu=save_lu(lu_pos)
            ipos=save_ipos(lu_pos)
            textline=save_textline(lu_pos)
            lu_pos=lu_pos-1
            write(*,*) 'returning to lu ',lu
          endif
#else    
          close(lu)
          lu=save_lu(lu_pos)
          ipos=save_ipos(lu_pos)
          textline=save_textline(lu_pos)
          lu_pos=lu_pos-1
          write(*,*) 'returning to lu ',lu
#endif
        else if (cval.eq.'quit'.or.cval.eq.'exit'.or.
     .           cval.eq.'end') then  ! quit/exit
          call terminate_plotting
          write(*,'(a,6i4)') 'STACK i,r,s,m,w,t = ',
     .     iistack,irstack,isstack,iwrk,iwwrk,itwrk
          call subend ()
          return
        else
          write(*,*) 'Unknown ',trim(cval) 
        endif  
      endif  
      goto 10

      contains

      subroutine select_matrix

      cflag=.true.
      CVALUE: select case (cval)
        case('pit','pitch')
          description='Pitch'
          nz=1
          do ix = -1, nx
            do iy = -1, ny
              if(bb(ix,iy,2).ne.0.0_R8) then
                mach(ix,iy,0)=
     .           sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/abs(bb(ix,iy,2))
              else
                mach(ix,iy,0)=1.0_R8
              endif
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case('bx','bpol')
          description='bx (T) (poloidal)'
          nz=1
          call copy(bb(:,:,0),1.0_R8,nz,textbl)
        case('by','brad')
          description='by (T)'
          nz=1
          call copy(bb(:,:,1),1.0_R8,nz,textbl)
        case('bz','btor')
          description='bz (T) (toroidal)'
          nz=1
          call copy(bb(:,:,2),1.0_R8,nz,textbl)
        case('bb')
          description='bb (T)'
          nz=1
          call copy(bb(:,:,3),1.0_R8,nz,textbl)
        case('b')
          description='B (T) (components)'
          nz=4
          call copy(bb(:,:,0),1.0_R8,nz,textcomponents)
        case('na')
          description='Ion density (m^{-3})'
          nz=ns
          call copy(na,1.0_R8,nz,textpl)
        case('nh')
          description='Hydrogenic atom density (m^{-3})'
          nz=1
          call copy(nh,1.0_R8,nz,textbl)
        case('ua')
          description='Parallel velocity (ua) (m.s^{-1})'
          nz=ns
          call copy(ua,1.0_R8,nz,textpl)
        case('ue')
          description='Parallel electron velocity (ue) (m.s^{-1})'
          nz=1
          call copy(ue,1.0_R8,nz,textbl)
        case('vv')
          description='Radial velocity (vv) (m.s^{-1})'
          nz=ns
          call copy(vv,1.0_R8,nz,textpl)
        case('ww')
          description='Toroidal velocity (ww) (m.s^{-1})'
          nz=ns
          do is=0,ns-1
            do ix=-1,nx
              do iy=-1,ny
                ww(ix,iy,is)=ua(ix,iy,is)*bb(ix,iy,2)/bb(ix,iy,3)-
     .                  uadia(ix,iy,0,is)*bb(ix,iy,0)/bb(ix,iy,3)
              end do
            end do
          end do
          call copy(ww,1.0_R8,nz,textpl)
        case('te')
          description='Electron temperature (eV)'
          nz=1
          call copy(te,1.0_R8/ev,nz,textbl)
        case('ti')
          description='Ion temperature (eV)'
          nz=1
          call copy(ti,1.0_R8/ev,nz,textbl)
        case('pr')
          description='Pressure (Pa)'
          nz=1
          call copy(pr,1.0_R8,nz,textbl)
        case('ne')
          description='Electron density (m^{-3})'
          nz=1
          call copy(ne,1.0_R8,nz,textbl)
        case('po')
          description='phi (V)'
          nz=1
          call copy(po,1.0_R8,nz,textbl)
        case('dfhe')
          description='diff.fhe (W)'
          nz=1
          call copy(div_fhe,1.0_R8,nz,textbl)
        case('dfhi')
          description='diff.fhi (W)'
          nz=1
          call copy(div_fhi,1.0_R8,nz,textbl)
        case('dhet','dfhet')
          description='diff.fhet (W)'
          nz=1
          call copy(div_fhet,1.0_R8,nz,textbl)
        case('dhit','dfhit')
          description='diff.fhit (W)'
          nz=1
          call copy(div_fhit,1.0_R8,nz,textbl)
        case('dfhj')
          description='diff.fhj (W)'
          nz=1
          call copy(div_fhj,1.0_R8,nz,textbl)
        case('dfch')
          description='diff.fch (A)'
          nz=1
          call copy(div_fch,1.0_R8,nz,textbl)
        case('dfna')
          description='div.fna (particles/s)'
          nz=ns
          call copy(div_fna,1.0_R8,nz,textpl)
        case('dfmo')
          description='div.fmo (N.s)'
          nz=ns
          call copy(div_fmo,1.0_R8,nz,textpl)
        case('dfhm')
          description='div.fhm (W)'
          nz=ns
          call copy(div_fhm,1.0_R8,nz,textpl)
        case('dfhp')
          description='div.fhp (W)'
          nz=ns
          call copy(div_fhp,1.0_R8,nz,textpl)
        case('zeff')
          description='Zeff'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              zeff(ix,iy)=0.0_R8
              do is=0,ns-1
                zeff(ix,iy)=zeff(ix,iy)+rz2(ix,iy,is)*na(ix,iy,is)
              end do
              zeff(ix,iy)=zeff(ix,iy)/ne(ix,iy)
            end do
          end do
          call copy(zeff,1.0_R8,nz,textbl)
        case('up')
          description='Poloidal velocity (up) (m.s^{-1})'
          nz=ns
          call copy(up,1.0_R8,nz,textpl)
        case('fnax')
          description='Poloidal particle flux (s^{-1})'
          nz=ns
          call copy(fnax,1.0_R8,nz,textpl)
        case('fnay')
          description='Radial particle flux (s^{-1})'
          nz=ns
          call copy(fnay,1.0_R8,nz,textpl)
        case('fmox')
          description='Poloidal flux of parallel momentum (N)'
          nz=ns
          call copy(fmox,1.0_R8,nz,textpl)
        case('fmoy')
          description='Radial flux of parallel momentum (N)'
          nz=ns
          call copy(fmoy,1.0_R8,nz,textpl)
        case('fhex')
          description='Poloidal electron energy flux (W)'
          nz=1
          call copy(fhe(:,:,0),1.0_R8,nz,textbl)
        case('fhey')
          description='Radial electron energy flux (W)'
          nz=1
          call copy(fhe(:,:,1),1.0_R8,nz,textbl)
        case('fhetx')
          description='Poloidal electron energy flux [fhetx] (W)'
          nz=1
          call copy(fhet(:,:,0),1.0_R8,nz,textbl)
        case('fhety')
        description='Radial electron energy flux [fhety] (W)'
          nz=1
          call copy(fhet(:,:,1),1.0_R8,nz,textbl)
        case('fhix')
          description='Poloidal ion energy flux (W)'
          nz=1
          call copy(fhi(:,:,0),1.0_R8,nz,textbl)
        case('fhiy')
          description='Radial ion energy flux (W)'
          nz=1
          call copy(fhi(:,:,1),1.0_R8,nz,textbl)
        case('fhitx')
          description='Poloidal ion energy flux [fhitx] (W)'
          nz=1
          call copy(fhit(:,:,0),1.0_R8,nz,textbl)
        case('fhity')
          description='Radial ion energy flux [fhity] (W)'
          nz=1
          call copy(fhit(:,:,1),1.0_R8,nz,textbl)
         case('fhmx')
          description='Poloidal parallel kinetic energy flux (W)'
          nz=ns
          call copy(fhmx,1.0_R8,nz,textpl)
        case('fhmy')
          description='Radial parallel kinetic energy flux (W)'
          nz=ns
          call copy(fhmy,1.0_R8,nz,textpl)
        case('fhpx')
          description='Poloidal potential energy flux (W)'
          nz=ns
          call copy(fhpx,1.0_R8,nz,textpl)
        case('fhpy')
          description='Radial potential energy flux (W)'
          nz=ns
          call copy(fhpy,1.0_R8,nz,textpl)
        case('fhjx')
          description='Poloidal electrostatic energy flux (W)'
          nz=1
          call copy(fhj(:,:,0),1.0_R8,nz,textbl)
        case('fhjy')
          description='Radial electrostatic energy flux (W)'
          nz=1
          call copy(fhj(:,:,1),1.0_R8,nz,textbl)
         case('fhtx')
          description='Poloidal total energy flux (W)'
          nz=1
          fhtx(:,:)=fhe(:,:,0)+fhi(:,:,0)+fhj(:,:,0)
          do is=0,ns-1
            fhtx(:,:)=fhtx(:,:)+fhp(:,:,0,is)
     &                         +fhm(:,:,0,is)*(1.0_R8-BoRiS)
          enddo
          call copy(fhtx,1.0_R8,nz,textbl)
        case('fhty')
          description='Radial total energy flux (W)'
          nz=1
          fhty(:,:)=fhe(:,:,1)+fhi(:,:,1)+fhj(:,:,1)
          do is=0,ns-1
            fhty(:,:)=fhty(:,:)+fhp(:,:,1,is)
     &                         +fhm(:,:,1,is)*(1.0_R8-BoRiS)
          enddo
          call copy(fhty,1.0_R8,nz,textbl)
        case('fle')
          description='Electron heat flux limit factor'
          nz=1
          call copy(fllime,1.0_R8,nz,textbl)
        case('fli')
          description='Ion heat flux limit factor'
          nz=1
          call copy(fllimi,1.0_R8,nz,textbl)
        case('nflfnax')
          description='fllim0fna_x'
          nz=ns
          call copy(fllim0fna_x,1.0_R8,nz,textpl)
        case('nflfnay')
          description='fllim0fna_y'
          nz=ns
          call copy(fllim0fna_y,1.0_R8,nz,textpl)
        case('nflfhix')
          description='fllim0fhi_x'
          nz=ns
          call copy(fllim0fhi_x,1.0_R8,nz,textpl)
        case('nflfhiy')
          description='fllim0fhi_y'
          nz=ns
          call copy(fllim0fhi_y,1.0_R8,nz,textpl)
        case('flvisc')
          description='fllimvisc'
          nz=ns
          call copy(fllimvisc,1.0_R8,nz,textpl)
        case('cvsx')
          description='cvsa_cl_x'
          nz=ns
          call copy(cvsa_cl_x,1.0_R8,nz,textpl)
        case('cvsy')
          description='cvsa_cl_y'
          nz=ns
          call copy(cvsa_cl_y,1.0_R8,nz,textpl)
        case('ceqp')
          description='ceqp'
          nz=1
	  write(*,*) 'ceqp min/max ',minval(ceqp),maxval(ceqp)
          call copy(ceqp,1.0_R8,nz,textbl)
       case('ebld')
          description='Energy balance (-divergence of fluxes) (W)'
          nz=1
          call copy(energy_balance(:,:),1.0_R8,nz,textbl)
       case('ebls')
          description='Energy balance (sources) (W)'
          nz=1
          call copy(energy_balance2(:,:),1.0_R8,nz,textbl)
       case('ebal')
          description='Energy balance (W)'
          nz=1
          call copy(energy_balance+energy_balance2,1.0_R8,nz,textbl)
       case('pbal')
          description='Particle balance (particles)'
          nz=ns
          call copy(particle_balance(:,:,:),1.0_R8,nz,textpl)
       case('hx')
          description='hx (m)'
          nz=1
          call copy(hx,1.0_R8,nz,textbl)
        case('hy')
          description='hy (m)'
          nz=1
          call copy(hy,1.0_R8,nz,textbl)
        case('hy1')
          description='hy1 (m) (hy corrected)'
          nz=1
          call copy(hy1,1.0_R8,nz,textbl)
        case('sx')
          description='sx (m^{2})'
          nz=1
          call copy(gs(:,:,0),1.0_R8,nz,textbl)
        case('sxperp')
          description='sxperp (m^{2})'
          nz=1
          call copy(gs(:,:,0)*qc(:,:),1.0_R8,nz,textbl)
        case('sy')
          description='sy (m^{2})'
          nz=1
          call copy(gs(:,:,1),1.0_R8,nz,textbl)
        case('vol')
          description='vol (m^{3})'
          nz=1
          call copy(vol,1.0_R8,nz,textbl)
        case('qz0')
          description='qz0'
          nz=1
          call copy(qz(:,:,0),1.0_R8,nz,textbl)
        case('qz1')
          description='qz1'
          nz=1
          call copy(qz(:,:,1),1.0_R8,nz,textbl)
        case('qc')
          description='qc'
          nz=1
          call copy(qc(:,:),1.0_R8,nz,textbl)
        case('mach')
          description='Mach number'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              vth=0.0_R8
              vel=0.0_R8
              do is=0,ns-1
                vth=vth+na(ix,iy,is)*
     .           (rza(ix,iy,is)*te(ix,iy)+ti(ix,iy))/
     .           (am(is)*ni(ix,iy,0)*mp)
                vel=vel+ua(ix,iy,is)*rza(ix,iy,is)*na(ix,iy,is)/
     .            ne(ix,iy)
              end do
              vth=sqrt(vth)
              mach(ix,iy,0)=vel/vth
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case('mach4')
          description='Mach number (SOLPS4.0)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mass=0.0_R8
              vth=0.0_R8
              vel=0.0_R8
              press=0.0_R8
              do is=0,ns-1
                if (.not.is_neutral(is)) then
                  mass=mass+am(is)*mp*na(ix,iy,is)
                  vel=vel+rza(ix,iy,is)*na(ix,iy,is)*ua(ix,iy,is)
                  press=press+na(ix,iy,is)*ti(ix,iy)
                endif
              end do
              press=press+ne(ix,iy)*te(ix,iy)
              vth=sqrt(press/mass)
              mach(ix,iy,0)=vel/(vth*ne(ix,iy))
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case('ehx')
          description='ehx'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ge.-1.and.rightix(ix,iy).le.nx.and.
     *         bottomiy(ix,iy).ge.-1.and.topiy(ix,iy).le.ny) then
                mach(ix,iy,0)=
     &           -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+
     &            (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))+
     &            (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))*0.5_R8*
     &            (log(ne(ix,iy))-log(ne(leftix(ix,iy),leftiy(ix,iy))))
     &            )/qe
              else
                mach(ix,iy,0)=0.0_R8
              endif
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case('ehy')
          description='ehy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ge.-1.and.rightix(ix,iy).le.nx.and.
     *         bottomiy(ix,iy).ge.-1.and.topiy(ix,iy).le.ny) then
                mach(ix,iy,0)=
     &           -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))+
     &            (ne(ix,iy)*te(ix,iy)-
     &             ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &             te(bottomix(ix,iy),bottomiy(ix,iy)))*
     &            (2.0_R8/qe)/
     &            (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))
              else
                mach(ix,iy,0)=0.0_R8
              endif
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case('dab2')
          description='Atomic density (m^{-3})'
          nz=natmi
          call copy(dab2(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textan)
        case('dmb2')
          description='Molecular density (m^{-3})'
          nz=nmoli
          call copy(dmb2(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textmn)
        case('ahal')
          description=
     1     'Atomic H-\Ga\R emissivity (photons.m^{-3}.s^{-1})'
          nz=1
          call copy(emiss(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textbl)
        case ('mhal')
          description=
     1     'Molecular H-\Ga\R emissivity (photons.m^{-3}.s^{-1})'
          nz=1
          call copy(emissmol(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textbl)
        case ('rqae')
          description='Electron cooling rate (W) [rqahe]'
          nz=ns
          call copy(rqahe,1.0_R8,nz,textpl)
        case ('b2ra')
          description='Line radiation rate (W)'
          nz=ns
          call copy(rqrad,1.0_R8,nz,textpl)
        case ('b2br')
          description='Bremsstrahlung radiation rate (W)'
          nz=ns
          call copy(rqbrm,1.0_R8,nz,textpl)
        case ('rsai')
          description='Ionisation ion energy losses (W) [rsahi]'
          nz=ns
          call copy(rsahi,1.0_R8,nz,textpl)
        case ('rsan')
          description='Ionisation particle losses (s^{-1}) [rsana]'
          nz=ns
          call copy(rsana,1.0_R8,nz,textpl)
        case ('rsam')
          description=
     1	   'Ionisation parallel momentum losses (m.s^{-2}) [rsamo]'
          nz=ns
          call copy(rsamo,1.0_R8,nz,textpl)
        case ('rran')
          description='Recombination particle losses (s^{-1}) [rrana]'
          nz=ns
          call copy(rrana,1.0_R8,nz,textpl)
        case ('rrai')
          description='Recombination ion energy losses (W) [rrahi]'
          nz=ns
          call copy(rrahi,1.0_R8,nz,textpl)
        case ('rram')
          description=
     1	   'Recombination parallel momentum losses (m.s^{-2}) [rramo]'
          nz=ns
          call copy(rramo,1.0_R8,nz,textpl)
        case ('rcxn')
          description='Particle neutral losses (s^{-1}) [rcxna]'
          nz=ns
          call copy(rcxna,1.0_R8,nz,textpl)
        case ('rcxi')
          description='Ion energy neutral losses (W) [rcxhi]'
          nz=ns
          call copy(rcxhi,1.0_R8,nz,textpl)
        case ('rcxm')
          description=
     1	   'Ion parallel momentum neutral losses (m.s^{-2}) [rcxmo]'
          nz=ns
          call copy(rcxmo,1.0_R8,nz,textpl)
        case ('fetx')
          description='Total poloidal internal energy flux (W)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              fetx(ix,iy)=fhe(ix,iy,0)+fhi(ix,iy,0)
              do is=0,ns-1
                fetx(ix,iy)=fetx(ix,iy)+
     .           rpt(ix,iy,is)*ev*fna(ix,iy,0,is)+
     1           fhm(ix,iy,0,is)*(1.0_R8-BoRiS)
              end do
            end do
          end do
          call copy(fetx,1.0_R8,nz,textbl)
        case ('fety')
          description='Total radial internal energy flux (W)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              fety(ix,iy)=fhe(ix,iy,1)+fhi(ix,iy,1)
              do is=0,ns-1
                fety(ix,iy)=fety(ix,iy)+
     1           rpt(ix,iy,is)*ev*fna(ix,iy,1,is)+
     1           fhm(ix,iy,1,is)*(1.0_R8-BoRiS)
              end do
            end do
          end do
          call copy(fety,1.0_R8,nz,textbl)
        case ('tab2')
          description='Neutral atom temperature (eV)'
          nz=natmi
          call copy(tab2(0:nx+1,0:ny+1,1:nz,1),1.0_R8/ev,nz,textan)
        case ('tmb2')
          description='Neutral molecule temperature (eV)'
          nz=nmoli
          call copy(tmb2(0:nx+1,0:ny+1,1:nz,1),1.0_R8/ev,nz,textmn)
        case ('rfla')
          description='Radial atomic flux (s^{-1}.m^{-2})'
          nz=natmi
          call copy(rfluxa(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textan)
        case ('rflm')
          description='Radial molecular flux (s^{-1}.m^{-2})'
          nz=nmoli
          call copy(rfluxm(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textmn)
        case ('refa')
          description='Radial atomic energy flux (W.m^{-2})'
          nz=natmi
          call copy(refluxa(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textan)
        case ('refm')
          description='Radial molecular energy flux (W.m^{-2})'
          nz=nmoli
          call copy(refluxm(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textmn)
        case ('pfla')
          description='Poloidal atomic flux (s^{-1}.m^{-2})'
          nz=natmi
          call copy(pfluxa(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textan)
        case ('pflm')
          description='Poloidal molecular flux (s^{-1}.m^{-2})'
          nz=nmoli
          call copy(pfluxm(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textmn)
        case ('pefa')
          description='Poloidal atomic energy flux (W.m^{-2})'
          nz=natmi
          call copy(pefluxa(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textan)
        case ('pefm')
          description=
     1             'Poloidal molecular energy flux (W.m^{-2})'
          nz=nmoli
          call copy(pefluxm(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textmn)
        case ('mch2')
          description='Normalised parallel velocity'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=ua(ix,iy,is)/
     1                      (ti(ix,iy)/(mp*am(is)))**0.5_R8
              end do
            end do
          end do
          call copy(mach,1.0_R8,nz,textpl)
        case ('ro')
          description='mass density (kg.m^{-3})'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=am(is)*mp*na(ix,iy,is)
              end do
            end do
          end do
          call copy(mach,1.0_R8,nz,textpl)
        case ('b2espcr')
          description='Eirene species index'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=real(b2espcr(is))
              end do
            end do
          end do
          call copy(mach,1.0_R8,nz,textpl)
        case ('mi')
          description='Ion mass (kg)'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=am(is)*mp
              end do
            end do
          end do
          call copy(mach,1.0_R8,nz,textpl)
        case ('zi','za')
          description='Ion charge'
          nz=ns
          call copy(rza,1.0_R8,nz,textpl)
        case ('za2')
          description='Squared ion charge'
          nz=ns
          call copy(rz2,1.0_R8,nz,textpl)
        case ('zn')
          description='Nuclear charge'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=zn(is)
              end do
            end do
          end do
          call copy(mach,1.0_R8,nz,textpl)
        case ('fmx')
          description='Poloidal mass flux (kg/s)'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=fna(ix,iy,0,is)*am(is)*mp
              end do
            end do
          end do
          call copy(mach,1.0_R8,nz,textpl)
        case ('fmy')
          description='Radial mass flux (kg/s)'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=fna(ix,iy,1,is)*am(is)*mp
              end do
            end do
          end do
          call copy(mach,1.0_R8,nz,textpl)
        case ('p*')
          description='Static pressure (Pa)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mass=0.0_R8
              vel=0.0_R8
              do is=0,ns-1
                mass=mass+mp*am(is)*na(ix,iy,is)
                vel=vel+rza(ix,iy,is)*na(ix,iy,is)*ua(ix,iy,is)
              end do
              vth=sqrt(pr(ix,iy)/mass)
              mach(ix,iy,0)=vel/(vth*ne(ix,iy))
              mach(ix,iy,0)=pr(ix,iy)*(1.0_R8+mach(ix,iy,0)**2)
            end do
          end do
          sumjxa=0.0_R8
          sumjxi=0.0_R8
          suml=0.0_R8
          sumr=0.0_R8
          iy=0
          ix=0
          do while(rightix(ix,iy).eq.ix+1 .and. ix.lt.nx)
            ix=ix+1
          end do
          ixl=ix+1
          ix=ix+1
          do while(rightix(ix,iy).eq.ix+1 .and. ix.lt.nx)
            ix=ix+1
          end do
          ixr=ix
          if(nnreg(0).eq.8) then
            ix=nx
            do while(leftix(ix,iy).eq.ix-1 .and. ix.gt.-1)
              ix=ix-1
            end do
            ixl=ix-1
            ix=ix-1
            do while(leftix(ix,iy).eq.ix-1 .and. ix.gt.-1)
              ix=ix-1
            end do
            ixr=ix
            ix=0
            do while(rightix(ix,iy).ne.nx+1 .and. ix.lt.nx)
              ix=ix+1
            end do
            ixr=ix
            do while(leftix(ix,iy).ne.-2 .and. ix.lt.nx)
              ix=ix+1
            end do
            ixl=ix
          endif
          ix=jxa
          do iy=jsep+1,ny
            sumjxa=sumjxa+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
          end do
          if(nnreg(0).eq.8) then
            ix=jxi
            do iy=jsep+1,ny
              sumjxi=sumjxi+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
            end do
          endif
          do iy=-1,ny
            ix=0
            suml(1)=suml(1)+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
            ix=nx-1
            sumr(1)=sumr(1)+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
            if(nnreg(0).eq.8) then
              ix=ixl+1
              suml(2)=suml(2)+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
              ix=ixr-1
              sumr(2)=sumr(2)+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
            endif
          end do
          write(*,*) 'Midplane pressure: ',sumjxa
          write(*,*) 'Inner target:      ',suml(1)
          write(*,*) 'Outer target:      ',sumr(1)
          write(*,*) 'Pressure drop (inner & outer): ',
     1             suml(1)/sumjxa,sumr(1)/sumjxa
          if(nnreg(0).eq.8) then
            write(*,*) 'Second X-point: '
            write(*,*) 'Midplane pressure: ',sumjxi
            write(*,*) 'Inner target:      ',suml(2)
            write(*,*) 'Outer target:      ',sumr(2)
            write(*,*) 'Pressure drop (inner & outer): ',
     1             suml(2)/sumjxi,sumr(2)/sumjxi
          endif
          call copy(mach,1.0_R8,nz,textbl)
        case ('conn')  
! New treatment using connectivity arrays
          description='Connection length (m)'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              t1=hx(ix,iy)/
     1         (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
              ix1=leftix(ix,iy)
              do while(ix1.ne.ix.and.ix1.ne.-2)
                t1=t1+hx(ix1,iy)/
     1           (sqrt(bb(ix1,iy,0)**2+bb(ix1,iy,1)**2)/bb(ix1,iy,3))
                ix1=leftix(ix1,iy)
              end do
              t2=hx(ix,iy)/
     1         (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
              ix1=rightix(ix,iy)
              do while(ix1.ne.ix.and.ix1.ne.nx+1)
                t2=t2+hx(ix1,iy)/
     1           (sqrt(bb(ix1,iy,0)**2+bb(ix1,iy,1)**2)/bb(ix1,iy,3))
                ix1=rightix(ix1,iy)
              end do
              mach(ix,iy,0)=min(t1,t2)
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('wpe')  
          description='omega_pe (rad/s)'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              mach(ix,iy,0)=sqrt(ne(ix,iy)*qe**2/me/eps0)
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('wce')  
          description='omega_ce (rad/s)'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              mach(ix,iy,0)=qe*bb(ix,iy,3)/me
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('wci')  
          description='omega_ci (rad/s)'
          nz=ns
          do is=0,ns-1
            do ix=-1,nx
              do iy=-1,ny
                mach(ix,iy,is)=rza(ix,iy,is)*qe*bb(ix,iy,3)/mp/am(is)
              end do
            end do
          enddo
          call copy(mach,1.0_R8,nz,textpl)
        case ('vthe')  
          description='v_th_e (m/s)'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              mach(ix,iy,0)=sqrt(te(ix,iy)/me)
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('vthi')  
          description='v_th_i (m/s)'
          nz=ns
          do is=0,ns-1
            do ix=-1,nx
              do iy=-1,ny
                mach(ix,iy,is)=sqrt(ti(ix,iy)/mp/am(is))
              end do
            end do
          enddo
          call copy(mach,1.0_R8,nz,textpl)
        case ('valfven')  
          description='valfven (m/s)'
          nz=ns
          do is=0,ns-1
            do ix=-1,nx
              do iy=-1,ny
                mach(ix,iy,is)=bb(ix,iy,3)/
     1           sqrt(mu0*na(ix,iy,is)*am(is)*mp)
              end do
            end do
          enddo
          call copy(mach,1.0_R8,nz,textpl)
        case ('cr')
          description='CR (m)'
          nz=1
          call copy(cr,1.0_R8,nz,textbl)
        case ('cz')
          description='CZ (m)'
          nz=1
          call copy(cz,1.0_R8,nz,textbl)
        case ('na0')
          description='na at previous time step (m^{-3})'
          nz=ns
          call copy(na0,1.0_R8,nz,textpl)
        case ('te0')
          description='te at previous time step (eV)'
          nz=1
          call copy(te0,1.0_R8/ev,nz,textbl)
        case ('ti0')
          description='ti at previous time step (eV)'
          nz=1
          call copy(ti0,1.0_R8/ev,nz,textbl)
        case ('ua0')
          description='ua at previous time step (m.s^{-1})'
          nz=ns
          call copy(ua0,1.0_R8,nz,textpl)
        case ('po0')
          description='po at previous time step (V)'
          nz=1
          call copy(po0,1.0_R8,nz,textbl)
        case ('smpr')
          description='Electrical field force'
          nz=ns
          write(*,*) 'Electrical field force : Not yet available!'
c          call copy(smpr,1.0_R8,nz,textpl)
        case ('smpt')
          description='Thermal force'
          nz=ns
          write(*,*) 'Thermal force : Not yet available!' 
c          call copy(smth,1.0_R8,nz,textpl)
        case ('smfr')
          description='Friction force'
          nz=ns
          write(*,*) 'Friction force : Not yet available!'
c          call copy(smfr,1.0_R8,nz,textpl)
        case ('pot')
        description='Potential energy (eV) (cum. ionisation energy)'
          nz=ns
          call copy(rpt,1.0_R8,nz,textpl)
        case ('poti')
          description='Ionisation energy (eV)'
          nz=ns
          call copy(rpi,1.0_R8,nz,textpl)
        case ('pot*')
        description='Ionisation energy (eV) [shifted down one]'
          nz=ns
          do is=0,ns-2
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=rpi(ix,iy,is+1)
              end do
            end do
          end do
          mach(:,:,ns-1)=0.0_R8
          call copy(mach,1.0_R8,nz,textpl)
        case ('smp')
          description='Self pressure force'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                ix1=leftix(ix,iy)
                if (ix1.ne.-2) then
                  mach(ix,iy,is)=-0.5_R8*
     1             abs(wbbl(ix,iy,0))/wbbl(ix,iy,3)*
     2             gs(ix,iy,0)*qc(ix,iy)*
     2             (na(ix1,iy,is)*ti(ix1,iy)-na(ix,iy,is)*ti(ix,iy))
                else
                  mach(ix,iy,is)=0.0_R8
                end if
              end do
            end do
          end do
          call copy(mach,1.0_R8,nz,textpl)
        case ('nreg')
          description='Region'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(istack(iistack) .eq. 0) then
                mach(ix,iy,0)=region(ix,iy,0)
              else if(istack(iistack) .gt. 0) then
                if(region(ix,iy,0) .eq. istack(iistack)) then
                  mach(ix,iy,0)=1.0_R8
                else
                  mach(ix,iy,0)=0.0_R8
                end if
              else
                if(region(ix,iy,0) .eq. -istack(iistack)) then
                  mach(ix,iy,0)=0.0_R8
                else
                  mach(ix,iy,0)=1.0_R8
                end if
              end if
            end do
          end do
          write(*,*) '"region" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('xreg')
          description='X-directed flux region'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(istack(iistack) .eq. 0) then
                mach(ix,iy,0)=region(ix,iy,1)
              else if(istack(iistack) .gt. 0) then
                if(region(ix,iy,1) .eq. istack(iistack)) then
                  mach(ix,iy,0)=1.0_R8
                else
                  mach(ix,iy,0)=0.0_R8
                end if
              else
                if(region(ix,iy,1) .eq. -istack(iistack)) then
                  mach(ix,iy,0)=0.0_R8
                else
                  mach(ix,iy,0)=1.0_R8
                end if
              end if
            end do
          end do
          write(*,*) 'X-directed flux "region" option with ',
     .     ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('yreg')
          description='Y-directed flux region'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(istack(iistack) .eq. 0) then
                mach(ix,iy,0)=region(ix,iy,2)
              else if(istack(iistack) .gt. 0) then
                if(region(ix,iy,2) .eq. istack(iistack)) then
                  mach(ix,iy,0)=1.0_R8
                else
                  mach(ix,iy,0)=0.0_R8
                end if
              else
                if(region(ix,iy,2) .eq. -istack(iistack)) then
                  mach(ix,iy,0)=0.0_R8
                else
                  mach(ix,iy,0)=1.0_R8
                end if
              end if
            end do
          end do
          write(*,*) 'Y-directed flux "region" option with ',
     .     ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('sna')
          description='sna'
          nz=ns
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          isel=ipop(istack,iistack)
          if(isel.lt.0.or.isel.gt.1) then
            write(*,*) 'selected value of ',isel,' out of range 0..1'
            return
          endif
          description=trim(description) // text0123(isel)
          do is=0,ns-1
            mach(-1:nx,-1:ny,is)=sna(-1:nx,-1:ny,isel,is)
          enddo
          write(*,*) '"sna" option with ',isel
          call copy(mach,1.0_R8,nz,textpl)
        case ('sne')
          description='sne'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          isel=ipop(istack,iistack)
          if(isel.lt.0.or.isel.gt.1) then
            write(*,*) 'selected value of ',isel,' out of range 0..1'
            return
          endif
          description=trim(description) // text0123(isel)
          write(*,*) '"sne" option with ',isel
          call copy(sne(:,:,isel),1.0_R8,nz,textbl)
        case ('smo')
          description='smo'
          nz=ns
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          isel=ipop(istack,iistack)
          if(isel.lt.0.or.isel.gt.3) then
            write(*,*) 'selected value of ',isel,' out of range 0..3'
            return
          endif
          description=trim(description) // text0123(isel)
          do is=0,ns-1
            mach(:,:,is)=smo(:,:,isel,is)
          enddo
          write(*,*) '"smo" option with ',isel
          call copy(mach,1.0_R8,nz,textpl)
        case ('smq')
          description='smq'
          nz=ns
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          isel=ipop(istack,iistack)
          if(isel.lt.0.or.isel.gt.3) then
            write(*,*) 'selected value of ',isel,' out of range 0..3'
            return
          endif
          description=trim(description) // text0123(isel)
          do is=0,ns-1
            mach(:,:,is)=smq(:,:,isel,is)
          enddo
          write(*,*) '"smq" option with ',isel
          call copy(mach,1.0_R8,nz,textpl)
        case ('she')
          description='she'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          isel=ipop(istack,iistack)
          if(isel.lt.0.or.isel.gt.3) then
            write(*,*) 'selected value of ',isel,' out of range 0..3'
            return
          endif
          description=trim(description) // text0123(isel)
          write(*,*) '"she" option with ',isel
          call copy(she(:,:,isel),1.0_R8,nz,textbl)
        case ('shi')
          description='shi'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          isel=ipop(istack,iistack)
          if(isel.lt.0.or.isel.gt.3) then
            write(*,*) 'selected value of ',isel,' out of range 0..3'
            return
          endif
          description=trim(description) // text0123(isel)
          write(*,*) '"shi" option with ',isel
          call copy(shi(:,:,isel),1.0_R8,nz,textbl)
        case ('sch')
          description='sch'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          isel=ipop(istack,iistack)
          if(isel.lt.0.or.isel.gt.3) then
            write(*,*) 'selected value of ',isel,' out of range 0..3'
            return
          endif
          description=trim(description) // text0123(isel)
          write(*,*) '"sch" option with ',isel
          call copy(sch(:,:,isel),1.0_R8,nz,textbl)
        case ('smaf')
          description='smaf'
          nz=ns
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          isel=ipop(istack,iistack)
          if(isel.lt.0.or.isel.gt.3) then
            write(*,*) 'selected value of ',isel,' out of range 0..3'
            return
          endif
          description=trim(description) // text0123(isel)
          do is=0,ns-1
            mach(:,:,is)=b2npmo_smaf(:,:,isel,is)
          enddo
          write(*,*) '"smaf" option with ',isel
          call copy(mach,1.0_R8,nz,textpl)
        case ('smag')
          description='smag'
          nz=ns
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          isel=ipop(istack,iistack)
          if(isel.lt.0.or.isel.gt.3) then
            write(*,*) 'selected value of ',isel,' out of range 0..3'
            return
          endif
          description=trim(description) // text0123(isel)
          do is=0,ns-1
            mach(:,:,is)=b2npmo_smag(:,:,isel,is)
          enddo
          write(*,*) '"smag" option with ',isel
          call copy(mach,1.0_R8,nz,textpl)
        case ('smav')
          description='smav'
          nz=ns
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          isel=ipop(istack,iistack)
          if(isel.lt.0.or.isel.gt.3) then
            write(*,*) 'selected value of ',isel,' out of range 0..3'
            return
          endif
          description=trim(description) // text0123(isel)
          do is=0,ns-1
            mach(:,:,is)=b2npmo_smav(:,:,isel,is)
          enddo
          write(*,*) '"smav" option with ',isel
          call copy(mach,1.0_R8,nz,textpl)
        case ('rco')
          description='RESCO'
          nz=ns
          call copy(resco,1.0_R8,nz,textpl)
        case ('rmo')
          description='RESMO'
          nz=ns
          call copy(resmo,1.0_R8,nz,textpl)
        case ('rmt')
          description='RESMT'
          nz=1
          call copy(resmt,1.0_R8,nz,textbl)
        case ('ree')
          description='RESEE'
          nz=1
          call copy(reshe,1.0_R8,nz,textbl)
        case ('rei')
          description='RESEI'
          nz=1
          call copy(reshi,1.0_R8,nz,textbl)
        case ('rpo')
          description='RESPO'
          nz=1
          call copy(respo,1.0_R8,nz,textbl)
        case ('cpa')
          description='CORPA'
          nz=ns
          call copy(corpa,1.0_R8,nz,textpl)
        case ('cua')
          description='CORUA'
          nz=ns
          call copy(corua,1.0_R8,nz,textpl)
        case ('cut')
          description='CORUT'
          nz=1
          call copy(corut,1.0_R8,nz,textbl)
        case ('cte')
          description='CORTE'
          nz=1
          call copy(corte,1.0_R8,nz,textbl)
        case ('cti')
          description='CORTI'
          nz=1
          call copy(corti,1.0_R8,nz,textbl)
        case ('cpo')
          description='CORPO'
          nz=1
          call copy(corpo,1.0_R8,nz,textbl)
        case ('ixlt')
          description='ixlt'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(ix .lt. istack(iistack)) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              end if
            end do
          end do
          write(*,*) '"ixlt" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('ixeq')
          description='ixeq'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(ix .eq. istack(iistack)) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              end if
            end do
          end do
          write(*,*) '"ixeq" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('ixgt')
          description='ixgt'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(ix .gt. istack(iistack)) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              end if
            end do
          end do
          write(*,*) '"ixgt" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('iylt')
          description='iylt'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(iy .lt. istack(iistack)) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              end if
            end do
          end do
          write(*,*) '"iylt" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('iyeq')
          description='iyeq'
          nz=1 
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(iy .eq. istack(iistack)) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              end if
            end do
          end do
          write(*,*) '"iyeq" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('iygt')
          description='iygt'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(iy .gt. istack(iistack)) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              end if
            end do
          end do
          write(*,*) '"iygt" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('islt')
          description='islt'
          nz=ns
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                if(is .lt. istack(iistack)) then
                  mach(ix,iy,is)=1.0_R8
                else
                  mach(ix,iy,is)=0.0_R8
                end if
              end do
            end do
          end do
          write(*,*) '"islt" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,ns,textpl)
        case ('iseq')
          description='iseq'
          nz=ns
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                if(is .eq. istack(iistack)) then
                  mach(ix,iy,is)=1.0_R8
                else
                  mach(ix,iy,is)=0.0_R8
                end if
              end do
            end do
          end do
          write(*,*) '"iseq" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,ns,textpl)
        case ('isgt')
          description='isgt'
          nz=ns
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                if(is .gt. istack(iistack)) then
                  mach(ix,iy,is)=1.0_R8
                else
                  mach(ix,iy,is)=0.0_R8
                end if
              end do
            end do
          end do
          write(*,*) '"isgt" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,ns,textpl)
        case ('lnlam')
          description='lnlam'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
!xpb lnlam is allowed to vary as per the NRL formulary using the
!xpb formula for ion-ion collisions or may be set by user.
!further modified to use Wesson
              if (lambda.eq.-5.0_R8) then                    !srv 09.02.11
                lnlam=max(-lambda,17.3_R8-
     1           log(ne(ix,iy)/1.0e20_R8)/2.0_R8+
     2           1.5_R8*log(ti(ix,iy)/ev/1000.0_R8))
              elseif (lambda.eq.-4.0_R8) then                !srv 09.02.11 {
	         if(te(ix,iy)/ev .le. 50.0_R8) then
                lnlam=23.4_R8 - 1.15_R8*log10(ne(ix,iy)/1.0e6_R8) +
     *                            3.45_R8*log10(te(ix,iy)/ev)
	         else
                lnlam=25.3_R8 - 1.15_R8*log10(ne(ix,iy)/1.0e6_R8) + 
     *                            2.30_R8*log10(te(ix,iy)/ev)
	         endif                                         !srv 09.02.11 }
              else
                lnlam = lambda
              endif
              mach(ix,iy,0)=lnlam
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('q')
          description='q'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              if(mod(region(ix,iy,0),4) .ne. 1) then
                mach(ix,iy,0)=0
              else
                mach(ix,iy,0)=safetyfactor(iy)
              end if
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('s')
          description='Shear'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              if(mod(region(ix,iy,0),4) .ne. 1) then
                mach(ix,iy,0)=0
              else
                mach(ix,iy,0)=shear(iy)
              end if
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('alph')
          description='\Ga\R'//
     1     ' (ratio or pressure gradient to critical ballooning value)'
          nz=1
          write(*,*) 'WARNING: "alph" not yet generalised'
c     assumes that bb is OK for 1d cases
          do iy=-1,ny
            do ix=-1,nx
              if(mod(region(ix,iy,0),4) .ne. 1) then
                mach(ix,iy,0)=0
              else
                if (iy.gt.-1) then
                  mach(ix,iy,0)=r0*
     1             (pr(ix,iy-1)-pr(ix,iy))/(rave(iy-1)-rave(iy))/
     2             (bb(ix,iy,3)**2/2.0_R8/mu0)/
     3             (0.0537_R8+0.534_R8*shear(iy)/safetyfactor(iy)**2)
                else
                  mach(ix,iy,0)=r0*
     1             (pr(ix,iy)-pr(ix,iy+1))/(rave(iy)-rave(iy+1))/
     2             (bb(ix,iy,3)**2/2.0_R8/mu0)/
     3             (0.0537_R8+0.534_R8*shear(iy)/safetyfactor(iy)**2)
                endif
              endif
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('d')
          description='D (m^{2}.s^{-1})'
          nz=ns
          call copy(dna0,1.0_R8,nz,textpl)
        case ('dp')
          description='Dpa (m^{2}.s^{-1})'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx     
                mach(ix,iy,is)=dpa0(ix,iy,is)*
     1           (rza(ix,iy,is)*te(ix,iy)+ti(ix,iy))
              enddo
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textpl)
        case ('kye')
          description='kye (m^{2}.s^{-1})'
          nz=1
          call copy(hce0(:,:)/ne(:,:),1.0_R8,nz,textbl)
        case ('kyi')
          description='kyi (m^{2}.s^{-1})'
          nz=ns
          call copy(hcib(:,:,:)/na(:,:,:),1.0_R8,nz,textpl)
        case ('kyi0')
          description='kyi0 (m^{2}.s^{-1})'
          nz=1
          call copy(hci0(:,:)/ni(:,:,0),1.0_R8,nz,textbl)
        case ('vlax')
          description='Poloidal pinch velocity (m.s^{-1})'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=vla0(ix,iy,0,is)
              enddo
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textpl)
        case ('vlay')
          description='Radial pinch velocity (m.s^{-1})'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=vla0(ix,iy,1,is)
              enddo
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textpl)
        case ('visc')
          description='Viscosity (m.kg^{-1}.s^{-1})'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=vsa0(ix,iy,is)/(mp*am(is)*na(ix,iy,is))
              enddo
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textpl)
        case ('bvis')
          description='Bas'' viscosity piece (m.kg^{-1}.s^{-1})'
          nz=ns
          call copy(bas_vis_piece,1.0_R8,nz,textpl)
        case ('sig')
          description='Anomalous conductivity (mho.m^{-1})'
          nz=1
          call copy(sig0(:,:)/(qe*ne(:,:)),1.0_R8,nz,textbl)
        case ('alf')
          description='Thermo-electric coefficient (mho.m^{-1})'
          nz=1
          call copy(alf0(:,:)/(ne(:,:)*sqrt(qe/te(:,:))),
     1     1.0_R8,nz,textbl)
        case ('fchp')
          description='Parallel current (A)'
          nz=1
          call copy(fch_p,1.0_R8,nz,textbl)
        case ('fchx')
          description='Poloidal current (A)'
          nz=1
          call copy(fch(:,:,0),1.0_R8,nz,textbl)
        case ('fchy')
          description='Radial current (A)'
          nz=1
          call copy(fch(:,:,1),1.0_R8,nz,textbl)
        case ('ud')
          description='pol comp of the ' //
     *          'drift velocity [uadia.x] (m.s^{-1})'
          nz=ns
          call copy(uadia(:,:,0,:),1.0_R8,nz,textpl)
        case ('vd')
          description='rad comp of the ' //
     *          'drift velocity [uadia.y] (m.s^{-1})'
          nz=ns
          call copy(uadia(:,:,1,:),1.0_R8,nz,textpl)
        case ('vadia.x')
          description='pol comp of the divergent ' //
     *          'part of the diamagnetic velocity [vadia.x] (m.s^{-1})'
          nz=ns
          call copy(uadia(:,:,0,:)-vaecrb(:,:,0,:),1.0_R8,nz,textpl)
        case ('vadia.y')
          description='rad comp of the divergent ' //
     *          'part of the diamagnetic velocity [vadia.y] (m.s^{-1})'
          nz=ns
          call copy(uadia(:,:,1,:)-vaecrb(:,:,1,:),1.0_R8,nz,textpl)
        case ('wadia.x')
          description='pol comp of the total ' //
     *          'diamagnetic velocity [wadia.x] (m.s^{-1})'
          nz=ns
          call copy(wadia(:,:,0,:),1.0_R8,nz,textpl)
        case ('wadia.y')
          description='rad comp of the total ' //
     *          'diamagnetic velocity [wadia.y] (m.s^{-1})'
          nz=ns
          call copy(wadia(:,:,1,:),1.0_R8,nz,textpl)
        case ('vexb.x')
          description='pol comp of the ' //
     *          'ExB velocity [vaecrb.x] (m.s^{-1})'
          nz=ns
          call copy(vaecrb(:,:,0,:),1.0_R8,nz,textpl)
        case ('vexb.y')
          description='rad comp of the ' //
     *          'ExB velocity [vaecrb.y] (m.s^{-1})'
          nz=ns
          call copy(vaecrb(:,:,1,:),1.0_R8,nz,textpl)
        case ('cd')
          description='fchdia.x'
          nz=1
          call copy(fchdia(:,:,0),1.0_R8,nz,textbl)
        case ('dd')
          description='fchdia.y'
          nz=1
          call copy(fchdia(:,:,1),1.0_R8,nz,textbl)
        case ('niso')
          description='Isolating region'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=0.0_R8
            enddo
          enddo
          do iso=1,nniso
            do iy=nyiso1(iso),nyiso2(iso)
              do ix=nxiso1(iso),nxiso2(iso)
                if(mach(ix,iy,0) .gt. 0.0_R8) then
                  write(*,*) 'Overlapping isolating region in cell',
     2                          ix,iy,'for isolating regions',
     3                          nint(mach(ix,iy,0)),'and',iso
                else
                  mach(ix,iy,0)=iso
                endif
              enddo
            enddo
           enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('bc')
  	  description='Boundary conditions'
	  nz=1
	  do iy=-1,ny
	    do ix=-1,nx
	      mach(ix,iy,0)=0.0_R8
	    enddo
	  enddo
	  do ibc=1,nbc
	    if(bcchar(ibc).eq.'S'.or.bcchar(ibc).eq.'N') then
	      iy=BCPOS(ibc)
	      do ix=BCSTART(ibc),BCEND(ibc)
		if(mach(ix,iy,0) .gt. 0.0_R8) then
		  write(*,*) 'Overlapping boundary conditions in cell ',
     2		   ix,iy,' for boundary condition ',
     3		   nint(mach(ix,iy,0)),'and',ibc
		else
		  mach(ix,iy,0)=ibc
		endif
	      enddo
	    else if(bcchar(ibc).eq.'E'.or.bcchar(ibc).eq.'W') then
	      ix=BCPOS(ibc)
	      do iy=BCSTART(ibc),BCEND(ibc)
		if(mach(ix,iy,0) .gt. 0.0_R8) then
		  write(*,*) 'Overlapping boundary conditions in cell ',
     2		   ix,iy,' for boundary condition ',
     3		   nint(mach(ix,iy,0)),'and',ibc
		else
		  mach(ix,iy,0)=ibc
		endif
	      enddo
	    else
	      write(*,*) 'Ignored BC with BCCHAR = ',bcchar(ibc),
     1	       ' for BC ',ibc
	    endif
	  enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('ix')
          description='ix (normalised)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(ix+1)/real(nx+1)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('iy')
          description='iy (normalised)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(iy+1)/real(ny+1)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('iix')
          description='ix'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(ix)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('iiy')
          description='iy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(iy)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('lix')
          description='leftix'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(leftix(ix,iy))
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('liy')
          description='leftiy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(leftiy(ix,iy))
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('rix')
          description='rightix'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(rightix(ix,iy))
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('riy')
          description='rightiy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(rightiy(ix,iy))
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('bix')
          description='bottomix'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(bottomix(ix,iy))
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('biy')
          description='bottomiy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(bottomiy(ix,iy))
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('tix')
          description='topix'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(topix(ix,iy))
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('tiy')
          description='topiy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(topiy(ix,iy))
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('rlcl')
          description='real cells'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ge.-1.and.leftix(ix,iy).le.nx.and.
     1         rightix(ix,iy).ge.-1.and.rightix(ix,iy).le.nx.and.
     1         bottomiy(ix,iy).ge.-1.and.bottomiy(ix,iy).le.ny.and.
     1         topiy(ix,iy).ge.-1.and.topiy(ix,iy).le.ny.and.
     2         region(ix,iy,0).ne.0) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              endif
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('gdcl')
          description='guard cells'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).lt.-1.or.leftix(ix,iy).gt.nx.or.
     1         rightix(ix,iy).lt.-1.or.rightix(ix,iy).gt.nx.or.
     1         bottomiy(ix,iy).lt.-1.or.bottomiy(ix,iy).gt.ny.or.
     1         topiy(ix,iy).lt.-1.or.topiy(ix,iy).gt.ny) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              endif
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('dib2')
          description='Test ion density (m^{-3})'
          nz=nioni
          call copy(dib2(0:nx+1,0:ny+1,1:nz,1),1.0_R8,nz,textin)
        case ('tib2')
          description='Test ion temperature (eV)'
          nz=nioni
          call copy(tib2(0:nx+1,0:ny+1,1:nz,1),1.0_R8/ev,nz,textin)
        case ('brna')
          description='b2stbr_sna'
          nz=ns
          call copy(b2stbr_sna,1.0_R8,nz,textpl)
        case ('brne')
          description='b2stbr_sne'
          nz=1
          call copy(b2stbr_sne,1.0_R8,nz,textbl)
        case ('brmo')
          description='b2stbr_smo'
          nz=ns
          call copy(b2stbr_smo,1.0_R8,nz,textpl)
        case ('brhe')
          description='b2stbr_she'
          nz=1
          call copy(b2stbr_she,1.0_R8,nz,textbl)
        case ('brhi')
          description='b2stbr_shi'
          nz=1
          call copy(b2stbr_shi,1.0_R8,nz,textbl)
        case ('brch')
          description='b2stbr_sch'
          nz=1
          call copy(b2stbr_sch,1.0_R8,nz,textbl)
        case ('bcna')
          description='b2stbc_sna'
          nz=ns
          call copy(b2stbc_sna,1.0_R8,nz,textpl)
        case ('bcne')
          description='b2stbc_sne'
          nz=1
          call copy(b2stbc_sne,1.0_R8,nz,textbl)
        case ('bcmo')
          description='b2stbc_smo'
          nz=ns
          call copy(b2stbc_smo,1.0_R8,nz,textpl)
        case ('bche')
          description='b2stbc_she'
          nz=1
          call copy(b2stbc_she,1.0_R8,nz,textbl)
        case ('bchi')
          description='b2stbc_shi'
          nz=1
          call copy(b2stbc_shi,1.0_R8,nz,textbl)
        case ('bcch')
          description='b2stbc_sch'
          nz=1
          call copy(b2stbc_sch,1.0_R8,nz,textbl)
        case ('exna')
          description='ext_sna'
          nz=ns
          call copy(ext_sna,1.0_R8,nz,textpl)
        case ('exne')
          description='ext_sne'
          nz=1
          call copy(ext_sne,1.0_R8,nz,textbl)
        case ('exmo')
          description='ext_smo'
          nz=ns
          call copy(ext_smo,1.0_R8,nz,textpl)
        case ('exhe')
          description='ext_she'
          nz=1
          call copy(ext_she,1.0_R8,nz,textbl)
        case ('exhi')
          description='ext_shi'
          nz=1
          call copy(ext_shi,1.0_R8,nz,textbl)
        case ('exch')
          description='ext_sch'
          nz=1
          call copy(ext_sch,1.0_R8,nz,textbl)
        case ('dvue')
          description='b2sihs_divue (phm0)'
          nz=1
          call copy(b2sihs_divue,1.0_R8,nz,textbl)
        case ('dvua')
          description='b2sihs_divua (phm1)'
          nz=1
          call copy(b2sihs_divua,1.0_R8,nz,textbl)
        case ('exbe')
          description='b2sihs_exbe (phm5)'
          nz=1
          call copy(b2sihs_exbe,1.0_R8,nz,textbl)
        case ('exba')
          description='b2sihs_exba (phm6)'
          nz=1
          call copy(b2sihs_exba,1.0_R8,nz,textbl)
        case ('visa')
          description='b2sihs_visa (phm2)'
          nz=1
          call copy(b2sihs_visa,1.0_R8,nz,textbl)
        case ('joul')
          description='b2sihs_joule (phm3)'
          nz=1
          call copy(b2sihs_joule,1.0_R8,nz,textbl)
        case ('fraa')
          description='b2sihs_fraa (phm4)'
          nz=1
          call copy(b2sihs_fraa,1.0_R8,nz,textbl)
        case ('str')
          description='b2sihs_str (style)'
          nz=1
          call copy(b2sihs_str,1.0_R8,nz,textbl)
        case ('ke')
          description='Parallel kinetic energy (J)'
          nz=ns
          call copy(kinrgy,1.0_R8,nz,textpl)
        case ('facdrift')
          description='facdrift'
          nz=1
          call copy(facdrift,1.0_R8,nz,textbl)
        case ('fac_exb')
          description='fac_ExB'
          nz=1
          call copy(fac_ExB,1.0_R8,nz,textbl)
        case ('fac_vis')
          description='fac_vis'
          nz=1
          call copy(fac_vis,1.0_R8,nz,textbl)
        case ('xymap')
          description='xymap'
          nz=1
          mach(-1:nx,-1:ny,0)=real(xymap(-1:nx,-1:ny))
          call copy(mach,1.0_R8,nz,textbl)
        case ('fnax32')
          description='fnax_32 (s^{-1})'
          nz=ns
          call copy(fna_32(:,:,0,:),1.0_R8,nz,textpl)
        case ('fnay32')
          description='fnay_32 (s^{-1})'
          nz=ns
          call copy(fna_32(:,:,1,:),1.0_R8,nz,textpl)
        case ('fnax52')
          description='fnax_52 (s^{-1})'
          nz=ns
          call copy(fna_52(:,:,0,:),1.0_R8,nz,textpl)
        case ('fnay52')
          description='fnay_52 (s^{-1})'
          nz=ns
          call copy(fna_52(:,:,1,:),1.0_R8,nz,textpl)
        case ('fnaxnodr')
          description='fnax_nodrift (s^{-1})'
          nz=ns
          call copy(fna_nodrift(:,:,0,:),1.0_R8,nz,textpl)
        case ('fnaynodr')
          description='fnay_nodrift (s^{-1})'
          nz=ns
          call copy(fna_nodrift(:,:,1,:),1.0_R8,nz,textpl)
        case ('fnix32')
          description='fnix_32 (s^{-1})'
          nz=1
          call copy(fni_32(:,:,0),1.0_R8,nz,textbl)
        case ('fniy32')
          description='fniy_32 (s^{-1})'
          nz=1
          call copy(fni_32(:,:,1),1.0_R8,nz,textbl)
        case ('fnix52')
          description='fnix_52 (s^{-1})'
          nz=1
          call copy(fni_52(:,:,0),1.0_R8,nz,textbl)
        case ('fniy52')
          description='fniy_52 (s^{-1})'
          nz=1
          call copy(fni_52(:,:,1),1.0_R8,nz,textbl)
        case ('fnex32')
          description='fnex_32 (s^{-1})'
          nz=1
          call copy(fne_32(:,:,0),1.0_R8,nz,textbl)
        case ('fney32')
          description='fney_32 (s^{-1})'
          nz=1
          call copy(fne_32(:,:,1),1.0_R8,nz,textbl)
        case ('fnex52')
          description='fnex_52 (s^{-1})'
          nz=1
          call copy(fne_52(:,:,0),1.0_R8,nz,textbl)
        case ('fney52')
          description='fney_52 (s^{-1})'
          nz=1
          call copy(fne_52(:,:,1),1.0_R8,nz,textbl)
        case ('fchx32')
          description='fchx_32 (A)'
          nz=1
          call copy(fch_32(:,:,0),1.0_R8,nz,textbl)
        case ('fchy32')
          description='fchy_32 (A)'
          nz=1
          call copy(fch_32(:,:,1),1.0_R8,nz,textbl)
        case ('fchx52')
          description='fchx_52 (A)'
          nz=1
          call copy(fch_52(:,:,0),1.0_R8,nz,textbl)
        case ('fchy52')
          description='fchy_52 (A)'
          nz=1
          call copy(fch_52(:,:,1),1.0_R8,nz,textbl)
        case ('fchxvspa')
          description='Pol comp of fchvispar (A)'
          nz=1
          call copy(fchvispar(:,:,0),1.0_R8,nz,textbl)
        case ('fchyvspa')
          description='Rad comp of fchvispar (A)'
          nz=1
          call copy(fchvispar(:,:,1),1.0_R8,nz,textbl)
        case ('fchxvspe')
          description='Pol comp of fchvisper (A)'
          nz=1
          call copy(fchvisper(:,:,0),1.0_R8,nz,textbl)
        case ('fchyvspe')
          description='Rad comp of fchvisper (A)'
          nz=1
          call copy(fchvisper(:,:,1),1.0_R8,nz,textbl)
        case ('fchxn')
          description='Pol comp of fchin (A)'
          nz=1
          call copy(fchin(:,:,0),1.0_R8,nz,textbl)
        case ('fchyn')
          description='Rad comp of fchin (A)'
          nz=1
          call copy(fchin(:,:,1),1.0_R8,nz,textbl)
#ifdef NEO
*    neoclassical part                                               !cjk
        case ('dnneo')
          description='Neoclassical particle diffusivity Dn ' //
     &                '(m^{2}.s^{-1})'
          nz=ns+1
          call copy(dn_neo,1.0_R8,nz,textplel)
        case ('chineo')
          description='Neoclassical heat conductivity Chi ' //
     &                '(m^{2}.s^{-1})'
          nz=ns+1
          call copy(chi_neo,1.0_R8,nz,textplel)
        case ('vnneo')
          description='Neoclassical particle pinch vn ' //
     &                '(m.s^{-1})'
          nz=ns+1
          call copy(vn_neo,1.0_R8,nz,textplel)
        case ('vwneo')
          description='Neoclassical heat pinch vw ' //
     &                '(m.s^{-1})'
          nz=ns+1
          call copy(vw_neo,1.0_R8,nz,textplel)
        case ('lnn')
          description='Density gradient length (m)'
          nz=ns+1
          call copy(lnn,1.0_R8,nz,textplel)
        case ('ltt')
          description='Temperature gradient length (m)'
          nz=ns+1
          call copy(ltt,1.0_R8,nz,textplel)
        case ('nuab')
          description='Two species collision frequency (s^{-1})'
          nz=(ns+1)*(ns+1)
          call copy(nuab,1.0_R8,nz,textplel2)
        case ('nua')
          description='Total collision frequency (s^{-1})'
          nz=ns+1
          call copy(nua,1.0_R8,nz,textplel)
        case ('nustar')
          description='Normalized collisionality'
          nz=ns+1
          call copy(nustar,1.0_R8,nz,textplel)
        case ('vth')
          description='Thermal velocity (m.s^{-1})'
          nz=ns+1
          call copy(vtha,1.0_R8/1.e2,nz,textplel)
        case ('rhol')
          description='Full-B Larmor radius (m)'
          nz=ns+1
          call copy(rhol,1.0_R8,nz,textplel)
        case ('deltal')
          description='Ratio of Larmor radius to minimum ' //
     &                'gradient length'
          nz=ns+1
          call copy(deltal,1.0_R8,nz,textplel)
        case ('wb')
          description='Width of the banana orbit (m)'
          nz=ns+1
          call copy(wb,1.0_R8,nz,textplel)
        case ('deltaw')
          description='Ratio of banana width to minimum ' //
     &                'gradient length'
          nz=ns+1
          call copy(deltaw,1.0_R8,nz,textplel)
        case ('deltawn')
          description='Ratio of banana width to density '//
     &                'gradient length'
          nz=ns+1
          call copy(deltawn,1.0_R8,nz,textplel)
        case ('deltawt')
          description='Ratio of banana width to temperature ' //
     &                'gradient length'
          nz=ns+1
          call copy(deltawt,1.0_R8,nz,textplel)
        case ('omega')
          description='Gyro-frequency (rad.s^{-1})'
          nz=ns+1
          call copy(omega,1.0_R8,nz,textplel)
        case ('deltab')
          description='Ratio of collision frequency to ' //
     &                'gyro-frequency'
          nz=ns+1
          call copy(deltab,1.0_R8,nz,textplel)
        case ('deltapar')
          description='Ratio of perpendicular to parallel dynamics'
          nz=ns+1
          call copy(deltapar,1.0_R8,nz,textplel)
        case ('chhichii')  
          description='Chang-Hinton \Gv\R_{i} (m^{2}/s)'
          call chang_hinton(nxd,nyd,nsd,jsep,jxa,lambda,
     1     safetyfactor,cr,cz,crx,bb,ne,te,ti,na,mach)
          nz=ns
          call copy(mach,1.0_R8,nz,textpl)
*    end neoclassical part
#endif
        case ('crx')
          description='Corner radial coordinate (m)'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(istack(iistack).ge.0.and.istack(iistack).le.3) then
            call copy(crx(:,:,istack(iistack)),1.0_R8,nz,
     1       textcorner(istack(iistack)))
          else
            write (*,*) 'Invalid corner index value'
          endif
          call pop_stack(stack_req)
        case ('cry')
          description='Corner vertical coordinate (m)'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(istack(iistack).ge.0.and.istack(iistack).le.3) then
            call copy(cry(:,:,istack(iistack)),1.0_R8,nz,
     1       textcorner(istack(iistack)))
          else
            write (*,*) 'Invalid corner index value'
          endif
          call pop_stack(stack_req)
        case default CVALUE
          cflag=.false.
      END SELECT CVALUE

      if(cflag) then
        call customise_label(description,graphlabel(iwrk))
        if(b2plot_debug.ge.1) 
     >   write(*,*) 'Copied ',trim(cval),' to plot array'
      endif
      return

      end subroutine select_matrix

      subroutine select_stack_manipulation

      cflag=.true.
      FVALUE: select case (cval)
        case('sumx')    ! sum over x
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call sumx()
          graphlabel(iwrk)=graphlabel(iwrk)
     1                 (1:lnblnk(graphlabel(iwrk)))//' SUMX '
        case('sumy')    ! sum over y
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call sumy()
          graphlabel(iwrk)=graphlabel(iwrk)
     1                 (1:lnblnk(graphlabel(iwrk)))//' SUMY '
        case('f.x')     ! line plot versus x
          stack_req=(/1,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.ny) then
            call fx(istack(iistack))
          else
            write(*,*) 'Faulty iy = ',istack(iistack),' chosen for f.x!'
          endif
          call pop_stack(stack_req)
        case('f.y')     ! line plot versus y
          stack_req=(/1,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.nx) then
            call fy(istack(iistack))
          else
            write(*,*) 'Faulty ix = ',istack(iistack),' chosen for f.y!'
          endif
          call pop_stack(stack_req)
        case('intx')    ! integrate over x
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call intx()
        case('inty')    ! integrate over y
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call inty()
        case('f.lx')    ! line plot versus log x
          stack_req=(/1,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.ny) then
            call flx(istack(iistack))
          else
            write(*,*)'Faulty iy = ',istack(iistack),' chosen for f.lx!'
          endif
          call pop_stack(stack_req)
        case('f.ly')    ! line plot versus log y
          stack_req=(/1,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.nx) then
            call fly(istack(iistack))
          else
            write(*,*)'Faulty ix = ',istack(iistack),' chosen for f.ly!'
          endif
          call pop_stack(stack_req)
        case('cont')    ! plot contours
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call cont
        case('surf')    ! plot patch plot
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call surf(1)
        case('line')    ! 
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call line(istack(iistack))
          call pop_stack(stack_req)
        case('psx')     ! per sx
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call psx
        case('psxperp')     ! per sxperp
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call psxperp
        case('psy')     ! per sy
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call psy
        case('abs')     ! abs()
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1                          abs(work(ixx(ix,iy),iyy(ix,iy),is,iwrk))
              end do
            end do
          end do
          graphlabel(iwrk)=
     .      graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//' (ABS)'
        case('ddx')     ! d/dx
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                ixl=leftix(ix,iy)
                ixr=rightix(ix,iy)
                if(ixl.eq.-2) then
                  mach(ix,iy,is)=
     1             (work(ixr,iy,is,iwrk)-work(ix,iy,is,iwrk))/
     2             (0.5_R8*(hx(ix,iy)+hx(ixr,iy)))
                elseif(ixr.eq.ndim(1,iwrk)+1) then
                  mach(ix,iy,is)=
     1             (work(ix,iy,is,iwrk)-work(ixl,iy,is,iwrk))/
     2             (0.5_R8*(hx(ix,iy)+hx(ixl,iy)))
                else
                  mach(ix,iy,is)=
     1             (work(ixr,iy,is,iwrk)-work(ixl,iy,is,iwrk))/
     2             (hx(ix,iy)+0.5_R8*(hx(ixl,iy)+hx(ixr,iy)))
                endif
              end do
            end do
          end do
          work(:,:,:,iwrk)=mach
          graphlabel(iwrk)=
     .      graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//'/(m) (d/dx)'
        case('ddy')     ! d/dy
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                iyb=bottomiy(ix,iy)
                iyt=topiy(ix,iy)
                if(iyb.eq.-2) then
                  mach(ix,iy,is)=
     1             (work(ix,iyt,is,iwrk)-work(ix,iy,is,iwrk))/
     2             (0.5_R8*(hy1(ix,iy)+hy1(ix,iyt)))
                elseif(iyt.eq.ndim(2,iwrk)+1) then
                  mach(ix,iy,is)=
     1             (work(ix,iy,is,iwrk)-work(ix,iyb,is,iwrk))/
     2             (0.5_R8*(hy1(ix,iy)+hy1(ix,iyb)))
                else
                  mach(ix,iy,is)=
     1             (work(ix,iyt,is,iwrk)-work(ix,iyb,is,iwrk))/
     2             (hy1(ix,iy)+0.5_R8*(hy1(ix,iyb)+hy1(ix,iyt)))
                endif
              end do
            end do
          end do
          work(:,:,:,iwrk)=mach
          graphlabel(iwrk)=
     .      graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//'/(m) (d/dy)'
        case('sumz')    ! sum over ns
          stack_req=(/2,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          izstrt=istack(iistack-1)
          izend=istack(iistack)
          if(izstrt .eq. 0 .and. izend .eq. 0) then
            izstrt=0
            izend=ndim(3,iwrk)-1
          end if
          if(izend .eq. -1) izend=ndim(3,iwrk)-1
          if(izstrt .lt. 0 .or. izend .gt. ndim(3,iwrk)-1 .or.
     1       izstrt .gt. izend) then
            write(*,*) 'Parameters supplied to sumz were invalid - ',
     1                 'valid parameters are:'
            write(*,*) ' 0 <= ',izstrt,' <= ',izend,' <= ',
     1           ndim(3,iwrk)-1
            return
          end if
          do is=izstrt+1,izend
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),izstrt,iwrk)=
     1                       work(ixx(ix,iy),iyy(ix,iy),izstrt,iwrk)
     2                       +work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              end do
            end do
          end do
          do iy=-1,ndim(2,iwrk)
            do ix=-1,ndim(1,iwrk)
              work(ixx(ix,iy),iyy(ix,iy),0,iwrk)=
     1                    work(ixx(ix,iy),iyy(ix,iy),izstrt,iwrk)
            end do
          end do
          ndim(3,iwrk)=1
          if(izstrt .ne. izend) then
            graphlabel(iwrk)=
     1       graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //' SUM'//
     2       textwk(izstrt,iwrk)(1:lnblnk(textwk(izstrt,iwrk)))
     3       //' to'//
     4       textwk(izend,iwrk)(1:lnblnk(textwk(izend,iwrk)))
          else
            graphlabel(iwrk)=
     1       graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))) //
     2       textwk(izstrt,iwrk)(1:lnblnk(textwk(izstrt,iwrk)))
          end if
          textwk(0,iwrk)=textbl(0)
          call pop_stack(stack_req)
        case('pvol')    ! per vol
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call pvol
        case('pne')     ! per ne
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call pne
        case('m+')      ! matrix addition
          stack_req=(/0,0,0,2,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          if(ndim(3,iwrk-1) .ne. ndim(3,iwrk) .and.
     1       ndim(3,iwrk-1) .ne. 1 .and. ndim(3,iwrk) .ne. 1) then
            write(*,*) 'Mismatch in ndim(3)''s in operands'
            return
          end if
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)=
     1     work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk-1)-1),iwrk-1)+
     2     work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk)-1),iwrk)
              end do
            end do
          end do
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk-1)=
     1                       work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)
              end do
            end do
            if(ndim(3,iwrk-1) .gt. ndim(3,iwrk)) then
              textwk(is,iwrk-1)=textwk(is,iwrk-1)
            else
              textwk(is,iwrk-1)=textwk(is,iwrk)
            end if
          end do
          graphlabel(iwrk-1)='Special'
          ndim(3,iwrk-1)=max(ndim(3,iwrk-1),ndim(3,iwrk))
          iwrk=iwrk-1
        case('m-')      ! matrix subtraction
          stack_req=(/0,0,0,2,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          if(ndim(3,iwrk-1) .ne. ndim(3,iwrk) .and.
     1       ndim(3,iwrk-1) .ne. 1 .and. ndim(3,iwrk) .ne. 1) then
            write(*,*) 'Mismatch in ndim(3)''s in operands'
            return
          end if
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)=
     1     work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk-1)-1),iwrk-1)-
     2     work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk)-1),iwrk)
              end do
            end do
          end do
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk-1)=
     2                       work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)
              end do
            end do
            if(ndim(3,iwrk-1) .gt. ndim(3,iwrk)) then
              textwk(is,iwrk-1)=textwk(is,iwrk-1)
            else
              textwk(is,iwrk-1)=textwk(is,iwrk)
            end if
          end do
          graphlabel(iwrk-1)='Special'
          ndim(3,iwrk-1)=max(ndim(3,iwrk-1),ndim(3,iwrk))
          iwrk=iwrk-1
        case('m*')      ! matrix multiplication
          stack_req=(/0,0,0,2,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          if(ndim(3,iwrk-1) .ne. ndim(3,iwrk) .and.
     1       ndim(3,iwrk-1) .ne. 1 .and. ndim(3,iwrk) .ne. 1) then
            write(*,*) 'Mismatch in ndim(3)''s in operands'
            return
          end if
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)=
     1     work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk-1)-1),iwrk-1)*
     2     work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk)-1),iwrk)
              end do
            end do
          end do
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk-1)=
     1                       work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)
              end do
            end do
            if(ndim(3,iwrk-1) .gt. ndim(3,iwrk)) then
              textwk(is,iwrk-1)=textwk(is,iwrk-1)
            else
              textwk(is,iwrk-1)=textwk(is,iwrk)
            end if
          end do
          graphlabel(iwrk-1)='Special'
          ndim(3,iwrk-1)=max(ndim(3,iwrk-1),ndim(3,iwrk))
          iwrk=iwrk-1
        case('m/')      ! matrix division
          stack_req=(/0,0,0,2,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          if(ndim(3,iwrk-1) .ne. ndim(3,iwrk) .and.
     1       ndim(3,iwrk-1) .ne. 1 .and. ndim(3,iwrk) .ne. 1) then
            write(*,*) 'Mismatch in ndim(3)''s in operands'
            return
          end if
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),
     1                   min(is,ndim(3,iwrk)-1),iwrk) .ne. 0.0_R8) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)=
     1     work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk-1)-1),iwrk-1)/
     2     work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk)-1),iwrk)
                else
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)=0.0_R8
                end if
              end do
            end do
          end do
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk-1)=
     2                       work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)
              end do
            end do
            if(ndim(3,iwrk-1) .gt. ndim(3,iwrk)) then
              textwk(is,iwrk-1)=textwk(is,iwrk-1)
            else
              textwk(is,iwrk-1)=textwk(is,iwrk)
            end if
          end do
          graphlabel(iwrk-1)='Special'
          ndim(3,iwrk-1)=max(ndim(3,iwrk-1),ndim(3,iwrk))
          iwrk=iwrk-1
        case('rm+')     ! matrix + real
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1           work(ixx(ix,iy),iyy(ix,iy),is,iwrk)+rstack(irstack)
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('rm-')     ! matrix - real
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1           work(ixx(ix,iy),iyy(ix,iy),is,iwrk)-rstack(irstack)
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('rm*')     ! matrix * real
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1           work(ixx(ix,iy),iyy(ix,iy),is,iwrk)*rstack(irstack)
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('rm/')     ! matrix / real
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          if(rstack(irstack).eq.0.0) then
            write (*,*) 'Divide by zero, You lose'
            return
          endif
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1           work(ixx(ix,iy),iyy(ix,iy),is,iwrk)/rstack(irstack)
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('rm**')    ! matrix ** real
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                  if (rstack(irstack).ne.0.0_R8) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1             work(ixx(ix,iy),iyy(ix,iy),is,iwrk)**rstack(irstack)
                else
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=1.0_R8
                end if
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('drop')    ! drop top of matrix stack
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          iwrk=iwrk-1
        case('dup')     ! duplicate top of matrix stack
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          if(iwrk .ge. wklng) call dropstackbottom()
          iwrk=iwrk+1
          ndim(1,iwrk)=ndim(1,iwrk-1)
          ndim(2,iwrk)=ndim(2,iwrk-1)
          ndim(3,iwrk)=ndim(3,iwrk-1)
          do is=0,ndim(3,iwrk-1)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1                       work(ixx(ix,iy),iyy(ix,iy),is,iwrk-1)
              end do
            end do
            textwk(is,iwrk)=textwk(is,iwrk-1)
          end do
          graphlabel(iwrk)=graphlabel(iwrk-1)
        case('swap')    ! swap top two elements of matrix stack
          stack_req=(/0,0,0,2,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          call iswap(ndim(1,iwrk),ndim(1,iwrk-1))
          call iswap(ndim(2,iwrk),ndim(2,iwrk-1))
          call iswap(ndim(3,iwrk),ndim(3,iwrk-1))
          do is=0,ndim(3,iwrk-1)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                call rswap(work(ixx(ix,iy),iyy(ix,iy),is,iwrk),
     1                       work(ixx(ix,iy),iyy(ix,iy),is,iwrk-1))
              end do
            end do
            call strswap(textwk(is,iwrk),textwk(is,iwrk-1))
          end do
          call strswap(graphlabel(iwrk),graphlabel(iwrk-1))
        case('over')    ! m2 m1 --> m2 m1 m2
          stack_req=(/0,0,0,2,0,0/)
          if(bad_args(stack_req)) return ! 2900
          if(iwrk .ge. wklng) call dropstackbottom()
          iwrk=iwrk+1
          ndim(1,iwrk)=ndim(1,iwrk-2)
          ndim(2,iwrk)=ndim(2,iwrk-2)
          ndim(3,iwrk)=ndim(3,iwrk-2)
          do is=0,ndim(3,iwrk-2)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1           work(ixx(ix,iy),iyy(ix,iy),is,iwrk-2)
              end do
            end do
            textwk(is,iwrk)=textwk(is,iwrk-2)
          end do
          graphlabel(iwrk)=graphlabel(iwrk-2)
        case('list')    ! list contents of matrix stack
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do i = iwrk, 1, -1
            write(*,*) trim(graphlabel(i)), 
     .       (trim(textwk(is,i)),is=0,ndim(3,i)-1)
          end do
        case('xsel')    ! select column
          stack_req=(/1,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          ix1=istack(iistack)
          if(-1.gt.ix1.or.ix1.gt.ndim(1,iwrk)) then
            write(*,*) 'Faulty parameter ix = ',ix1,' for xsel!'
            return
          endif
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1           work(ixx(ix1,iy),iyy(ix1,iy),is,iwrk)
              end do
            end do
          end do
          call pop_stack(stack_req)
        case('ysel')    ! select row
          stack_req=(/1,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          iy1=istack(iistack)
          if(-1.gt.iy1.or.iy1.gt.ndim(2,iwrk)) then
            write(*,*) 'Faulty parameter iy = ',iy1,' for ysel!'
            return
          endif
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1           work(ixx(ix,iy1),iyy(ix,iy1),is,iwrk)
              end do
            end do
          end do
          call pop_stack(stack_req)
        case('zsel')    ! select species
          stack_req=(/1,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          is1=istack(iistack)
          if(0.gt.is1.or.is1.ge.ndim(3,iwrk)) then
            write(*,*) 'Faulty parameter is = ',is1,' for zsel!'
            return
          endif
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1                      work(ixx(ix,iy),iyy(ix,iy),is1,iwrk)
              end do
            end do
            textwk(is,iwrk)=textwk(is1,iwrk)
          end do
          ndim(3,iwrk)=1
          call pop_stack(stack_req)
        case('mmin')    ! min(m1,m2)
          stack_req=(/0,0,0,2,0,0/) 
          if(bad_args(stack_req)) return ! 2900
          if(ndim(3,iwrk-1) .ne. ndim(3,iwrk) .and.
     1       ndim(3,iwrk-1) .ne. 1 .and. ndim(3,iwrk) .ne. 1) then
            write(*,*) 'Mismatch in ndim(3)''s in operands'
            return
          end if
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)=
     1  min(work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk-1)-1),iwrk-1),
     2      work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk)-1),iwrk))
              end do
            end do
          end do
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk-1)=
     1                       work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)
              end do
            end do
            if(ndim(3,iwrk-1) .gt. ndim(3,iwrk)) then
              textwk(is,iwrk-1)=textwk(is,iwrk-1)
            else
              textwk(is,iwrk-1)=textwk(is,iwrk)
            end if
          end do
          graphlabel(iwrk-1)='Special'
          ndim(3,iwrk-1)=max(ndim(3,iwrk-1),ndim(3,iwrk))
          iwrk=iwrk-1
        case('mmax')    ! max(m1,m2)
          stack_req=(/0,0,0,2,0,0/)
          if(bad_args(stack_req)) return ! 2900
          if(ndim(3,iwrk-1) .ne. ndim(3,iwrk) .and.
     1       ndim(3,iwrk-1) .ne. 1 .and. ndim(3,iwrk) .ne. 1) then
            write(*,*) 'Mismatch in ndim(3)''s in operands'
            return
          end if
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)=
     1  max(work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk-1)-1),iwrk-1),
     2      work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk)-1),iwrk))
              end do
            end do
          end do
          do is=0,max(ndim(3,iwrk-1),ndim(3,iwrk))-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk-1)=
     1                       work(ixx(ix,iy),iyy(ix,iy),is,iwrk+1)
              end do
            end do
            if(ndim(3,iwrk-1) .gt. ndim(3,iwrk)) then
              textwk(is,iwrk-1)=textwk(is,iwrk-1)
            else
              textwk(is,iwrk-1)=textwk(is,iwrk)
            end if
          end do
          graphlabel(iwrk-1)='Special'
          ndim(3,iwrk-1)=max(ndim(3,iwrk-1),ndim(3,iwrk))
          iwrk=iwrk-1
        case('mlt')    ! m<r
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),is,iwrk).lt.
     1           rstack(irstack)) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=1.0_R8
                else
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=0.0_R8
                endif
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('mle')    ! m<=r
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),is,iwrk).le.
     1           rstack(irstack)) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=1.0_R8
                else
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=0.0_R8
                endif
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('meq')    ! m=r
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),is,iwrk).eq.
     1           rstack(irstack)) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=1.0_R8
                else
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=0.0_R8
                endif
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('mne')    ! m!=r
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),is,iwrk).ne.
     1           rstack(irstack)) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=1.0_R8
                else
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=0.0_R8
                endif
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('mgt')    ! m>r
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),is,iwrk).gt.
     1           rstack(irstack)) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=1.0_R8
                else
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=0.0_R8
                endif
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('mge')    ! m>=r
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),is,iwrk).ge.
     1           rstack(irstack)) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=1.0_R8
                else
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=0.0_R8
                endif
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('rmmi')    ! min(r1,m1)
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1        min(work(ixx(ix,iy),iyy(ix,iy),is,iwrk),rstack(irstack)) 
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('rmma')    ! max(r1,m1)
          stack_req=(/0,1,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1        max(work(ixx(ix,iy),iyy(ix,iy),is,iwrk),rstack(irstack))
              end do
            end do
          end do
          dummy=rpop(rstack,irstack)
          graphlabel(iwrk)='Special'
        case('rtom')    ! real --> matrix stack
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(iwrk .eq. wklng) call dropstackbottom()
          iwrk=iwrk+1
          ndim(1,iwrk)=nx
          ndim(2,iwrk)=ny
          ndim(3,iwrk)=1
          work(-1:nx,-1:ny,0,iwrk)=rstack(irstack)
          values(9)=rstack(irstack)
          values(10)=rstack(irstack)
          textwk(0,iwrk)=textbl(0)
          write(graphlabel(iwrk),'(1pe15.8)') rpop(rstack,irstack)
        case('log')     ! log10()
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),is,iwrk).gt.0.0_R8) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1                       log10(work(ixx(ix,iy),iyy(ix,iy),is,iwrk))
                else
                  write(*,*) 'Negative argument for ',
     1              ixx(ix,iy),iyy(ix,iy),is,iwrk,' : Nothing done.'
                end if
              end do
            end do
          end do
          graphlabel(iwrk)='Special'
        case('ln')      ! ln()
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),is,iwrk).gt.0.0_R8) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1                          log(work(ixx(ix,iy),iyy(ix,iy),is,iwrk))
                else
                  write(*,*) 'Negative argument for ',
     1        ixx(ix,iy),iyy(ix,iy),is,iwrk,'. Nothing done.'
                end if
              end do
            end do
          end do
          graphlabel(iwrk)='Special'
        case('sqrt')    ! sqrt()
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),is,iwrk).ge.0.0_R8) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1                       sqrt(work(ixx(ix,iy),iyy(ix,iy),is,iwrk))
                else
                  write(*,*) 'Negative argument for ',
     1              ixx(ix,iy),iyy(ix,iy),is,iwrk,' : Nothing done.'
                end if
              end do
            end do
          end do
          call pop_stack(stack_req)
          graphlabel(iwrk)='Special'
        case('exp')     ! exp()
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1                       exp(work(ixx(ix,iy),iyy(ix,iy),is,iwrk))
              end do
            end do
          end do
          graphlabel(iwrk)='Special'
        case('zrng')    ! select species range
          stack_req=(/2,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          if(0.gt.istack(iistack-1).or.
     &       istack(iistack-1).gt.ndim(3,iwrk)-1) then
            write(*,*) 'Faulty first parameter is = ',istack(iistack-1),
     &       ' for zrng!'
            return
          endif
          if(0.gt.istack(iistack).or.
     &       istack(iistack).gt.ndim(3,iwrk)-1) then
            write(*,*) 'Faulty second parameter is = ',istack(iistack),
     &       ' for zrng!'
            return
          endif
          do is=istack(iistack-1),istack(iistack)
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is-istack(iistack-1),iwrk)=
     1                       work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              end do
            end do
            textwk(is-istack(iistack-1),iwrk)=textwk(is,iwrk)
          end do
          ndim(3,iwrk)=istack(iistack)-istack(iistack-1)+1
          call pop_stack(stack_req)
        case('chs')     ! change sign
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1                       -work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              end do
            end do
          end do
          graphlabel(iwrk)=
     1      graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//' (CHS)'
        case('bala')    ! balance ???
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          if(.not.(ndim(3,iwrk) .eq. 1 .or. ndim(3,iwrk) .eq. ns)) then
            write(*,*) 'Balancing only works for dimensions 1 or ns'
            write(*,*) 'Setting dimension to 1'
            ndim(3,iwrk)=1
          endif
          if(ndim(3,iwrk) .eq. 1) then
            do is=1,ns-1
              do iy=-1,ndim(2,iwrk)
                do ix=-1,ndim(1,iwrk)
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1              work(ixx(ix,iy),iyy(ix,iy),0,iwrk)
                enddo
              enddo
            enddo
          endif
          do isp=0,nspecies-1
            balance(isp)=0.0_R8
          enddo
          do is=0,ns-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                balance(ispecies(is))=balance(ispecies(is))+
     1            work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              enddo
            enddo
          enddo
          do isp=0,nspecies-1
            balance(isp)=1.0_R8/balance(isp)
          enddo
          do is=0,ns-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1      work(ixx(ix,iy),iyy(ix,iy),is,iwrk)*balance(ispecies(is))
              enddo
            enddo
          enddo
        case('traj')    ! trajectories
          stack_req=(/5,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(istack(5).gt.0) then
            stack_req=(/7+istack(iistack),0,0,0,0,0/)
          else
            stack_req=(/7+istack(iistack),0,0,1,0,0/)
          endif
          if(bad_args(stack_req)) return
          call caltraj(istack(iistack-6-istack(iistack)),
     +     7+istack(iistack))
          call pop_stack(stack_req)
        case('stow')    ! matrix --> wall array
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          if(iwwrk .eq. wklng) call dropwallstackbottom()
          iwwrk=iwwrk+1
          nwdim(1,iwwrk)=nbnd
          nwdim(2,iwwrk)=ndim(3,iwrk)
          do i = 1, nbnd
            wall_stack(i,0:ndim(3,iwrk)-1,iwwrk)=
     .       work(imapx(i),imapy(i),0:ndim(3,iwrk)-1,iwrk)
          enddo
          values(14)=minval(wall_stack(1:nbnd,0:nwdim(2,iwwrk)-1,iwwrk))
          values(15)=maxval(wall_stack(1:nbnd,0:nwdim(2,iwwrk)-1,iwwrk))
          wall_label(iwwrk)=graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))
          textwwk(0:nwdim(2,iwwrk),iwwrk)=textwk(0:ndim(3,iwrk),iwrk)
          goto 2059
        case('xtow')    ! x dimension of matrix --> wall array
          stack_req=(/1,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          iy = ipop(istack,iistack)
          if(iwwrk .eq. wklng) call dropwallstackbottom()
          iwwrk=iwwrk+1
          nwdim(1,iwwrk)=ndim(1,iwrk)+2
          nwdim(2,iwwrk)=ndim(3,iwrk)
          do ix=-1,ndim(1,iwrk)
            wall_stack(ix+2,0:ndim(3,iwrk)-1,iwwrk)=
     .       work(ixx(ix,iy),iyy(ix,iy),0:ndim(3,iwrk)-1,iwrk)
          enddo
          values(14)=
     .     minval(wall_stack(1:nwdim(1,iwwrk),0:nwdim(2,iwwrk)-1,iwwrk))
          values(15)=
     .     maxval(wall_stack(1:nwdim(1,iwwrk),0:nwdim(2,iwwrk)-1,iwwrk))
          wall_label(iwwrk)=graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))
          textwwk(0:nwdim(2,iwwrk),iwwrk)=textwk(0:ndim(3,iwrk),iwrk)
          goto 2059
        case('ytow')    ! y dimension of matrix --> wall array
          stack_req=(/1,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          if(iwwrk .eq. wklng) call dropwallstackbottom()
          ix = ipop(istack,iistack)
          iwwrk=iwwrk+1
          nwdim(1,iwwrk)=ndim(2,iwrk)+2
          nwdim(2,iwwrk)=ndim(3,iwrk)
          do iy=-1,ndim(2,iwrk)
            wall_stack(iy+2,0:ndim(3,iwrk)-1,iwwrk)=
     .       work(ixx(ix,iy),iyy(ix,iy),0:ndim(3,iwrk)-1,iwrk)
          enddo
          values(14)=
     .     minval(wall_stack(1:nwdim(1,iwwrk),0:nwdim(2,iwwrk)-1,iwwrk))
          values(15)=
     .     maxval(wall_stack(1:nwdim(1,iwwrk),0:nwdim(2,iwwrk)-1,iwwrk))
          wall_label(iwwrk)=graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))
          textwwk(0:nwdim(2,iwwrk),iwwrk)=textwk(0:ndim(3,iwrk),iwrk)
          goto 2059
        case('rtow')    ! real --> wall stack
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(iwwrk .eq. wklng) call dropwallstackbottom()
          iwwrk=iwwrk+1
          nwdim(1,iwwrk)=nbnd
          nwdim(2,iwwrk)=1
          wall_stack(1:nbnd,0,iwwrk)=rstack(irstack)
          values(14)=rstack(irstack)
          values(15)=rstack(irstack)
          textwwk(0,iwwrk)=textbl(0)
          write(wall_label(iwwrk),'(1pe15.8)') rpop(rstack,irstack)
          call set_wallbarlabel(wall_label(iwwrk))
        case('rtot')    ! real --> target stack
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(itwrk .eq. wklng) call droptargetstackbottom()
          itwrk=itwrk+1
          target_stack(1:nbnd,1:ndepth_nml,itwrk)=rstack(irstack)
          values(32)=rstack(irstack)
          values(33)=rstack(irstack)
          write(target_label(itwrk),'(1pe15.8)') rpop(rstack,irstack)
          call set_targetbarlabel(target_label(itwrk))
        case('eirc')    ! add the Eirene atomic densities to the b2 values
          stack_req=(/0,0,0,0,0,0/)
          if(b2plot_debug.ge.12) 
     &     write(*,*) iistack,irstack,isstack,iwrk,iwwrk,stack_req
          if(bad_args(stack_req)) write(*,*) 'Sorry!'
          do is=1,natmi
             do iss=0,ns-1
                stmp1=trim(textpl(iss))
                call lower_case(stmp1)                                !srv 14.10.13
                stmp2=trim(textan(is-1))//'0'
                call lower_case(stmp2)                                !srv 14.10.13
                if(is_neutral(iss) .and. stmp1.eq.stmp2) then
                   write(*,*) 'We seem to have found a match between ',
     &                  trim(textpl(iss)),' and ',
     &                  trim(textan(is-1))//'0',
     &                  '  (',iss,' & ',is,')'
                   do iy=-1,ny
                      do ix=-1,nx
                         na(ix,iy,iss)=na(ix,iy,iss)+
     &                        dab2(ix+1,iy+1,is,1)
                      enddo
                   enddo
                endif
             enddo
          enddo
          return
        case('nese')    ! set ne from matrix stack
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do iy=-1,ndim(2,iwrk)
            do ix=-1,ndim(1,iwrk)
              ne(ix,iy)=work(ixx(ix,iy),iyy(ix,iy),0,iwrk)
            enddo
          enddo
          goto 2059
        case('tese')    ! set te from matrix stack
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ndim(2,iwrk)
            do ix=-1,ndim(1,iwrk)
              te(ix,iy)=work(ixx(ix,iy),iyy(ix,iy),0,iwrk)
            enddo
          enddo
          goto 2059        ! drop
        case('dase')    ! set dab2 from matrix stack
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,natmi-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                dab2(ix,iy,is,1)=
     .           work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk)-1),iwrk)
              enddo
            enddo
          enddo
          goto 2059
        case('uuse')    ! set uutraj from matrix stack
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ns-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                uutraj(ix,iy,is)=work(ixx(ix,iy),iyy(ix,iy),
     1                             min(is,ndim(3,iwrk)-1),iwrk)
              enddo
            enddo
          enddo
          goto 2059        ! drop
        case('vvse')    ! set vvtraj from matrix stack
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ns-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                vvtraj(ix,iy,is)=work(ixx(ix,iy),iyy(ix,iy),
     1                            min(is,ndim(3,iwrk)-1),iwrk)
              enddo
            enddo
          enddo
          goto 2059
        case('arox')    ! set arrow field x-component from matrix
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ns-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                arrowx(ix,iy,is)=
     1           work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk)-1),iwrk)
              enddo
            enddo
          enddo
          goto 2059
        case('aroy')    ! set arrow field y-component from matrix
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do is=0,ns-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                arrowy(ix,iy,is)=
     1           work(ixx(ix,iy),iyy(ix,iy),min(is,ndim(3,iwrk)-1),iwrk)
              enddo
            enddo
          enddo
          goto 2059
        case('pos')     ! matrix with 1 if element > 0, else 0
          stack_req=(/0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                if(work(ixx(ix,iy),iyy(ix,iy),is,iwrk) .gt. 0.0_R8) then
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=1.0_R8
                else
                  work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=0.0_R8
                endif
              end do
            end do
          end do
          graphlabel(iwrk)='Special'
        case('target')   ! add target elements plot to surface plot
          call targeton(.true.,values(32),values(33),trglinlog)
        case('trgoff')   ! turns off target elements plot
          call targetoff 
        case('idrop')   ! drop top of integer stack
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp=ipop(istack,iistack)
        case('idup')    ! duplicate top of integer stack
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp=istack(iistack)
          call ipush(istack,iistack,itmp)
        case('iswap')   ! swap top two elements of integer stack
          stack_req=(/2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call iswap(istack(iistack),istack(iistack-1))
        case('iover')   ! i2 i1 --> i2 i1 i2
          stack_req=(/2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp=istack(iistack-1)
          call ipush(istack,iistack,itmp)
        case('rdrop')   ! drop top of real stack
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          t1=rpop(rstack,irstack)
        case('rdup')    ! duplicate top of real stack
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          t1=rstack(irstack)
          call rpush(rstack,irstack,t1)
        case('rswap')   ! swap top two elements of real stack
          stack_req=(/0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          call rswap(rstack(irstack),rstack(irstack-1))
        case('rover')   ! r2 r1 --> r2 r1 r2
          stack_req=(/0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          t1=rstack(irstack-1)
          call rpush(rstack,irstack,t1)
        case default FVALUE
          cflag=.false.
      END SELECT FVALUE
c
      goto 2900

 2059 continue
      stack_req=(/0,0,0,1,0,0/)
      if(bad_args(stack_req)) return
      iwrk=iwrk-1

 2900 continue
      if(cflag) then
        if(b2plot_debug.ge.4)
     1   write(*,*) 'Applied operator ',trim(cval),' to plot array'
        if(iwrk.gt.0) then
          values(9)=minval(work(-1:ndim(1,iwrk),-1:ndim(2,iwrk),
     1         0:ndim(3,iwrk)-1,iwrk))
          values(10)=maxval(work(-1:ndim(1,iwrk),-1:ndim(2,iwrk),
     1         0:ndim(3,iwrk)-1,iwrk))
        endif
      endif
      return

      end subroutine select_stack_manipulation

      subroutine select_parameters

      cflag=.true.
      stack_req=(/0,1,0,0,0,0/)
      NVALUE: select case(cval)
        case('xmin')    ! Minimum x for plotting in the computational domain
          if(bad_args(stack_req)) return
          values(1)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(1)
        case('xmax')    ! Maximum x for plotting in the computational domain
          if(bad_args(stack_req)) return
          values(2)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(2)
        case('ymin')    ! Minimum y for plotting in the computational domain
          if(bad_args(stack_req)) return
          values(3)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(3)
        case('ymax')    ! Maximum y for plotting in the computational domain
          if(bad_args(stack_req)) return
          values(4)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(4)
        case('rmin')    ! Minimum r for plotting in the physical domain
          if(bad_args(stack_req)) return
          values(5)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(5)
        case('rmax')    ! Maximum r for plotting in the physical domain
          if(bad_args(stack_req)) return
          values(6)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(6)
        case('zmin')    ! Minimum z for plotting in the physical domain
          if(bad_args(stack_req)) return
          values(7)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(7)
        case('zmax')    ! Maximum z for plotting in the physical domain
          if(bad_args(stack_req)) return
          values(8)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(8)
        case('fmin')    ! Minimum function value for plotting
          if(bad_args(stack_req)) return
          values(9)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(9)
        case('fmax')    ! Maximum function value for plotting
          if(bad_args(stack_req)) return
          values(10)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(10)
        case('omin')    ! Minimum ordinate for 2-d plots
          if(bad_args(stack_req)) return
          values(11)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(11)
        case('omax')    ! Maximum ordinate for 2-d plots
          if(bad_args(stack_req)) return
          values(12)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(12)
        case('xref')    ! Reference position in x: 0.0=left, 0.5=mid, 1.0=right
          if(bad_args(stack_req)) return
          values(13)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(13)
        case('wmin')    ! wmin in pltwall
          if(bad_args(stack_req)) return
          values(14)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(14)
        case('wmax')    ! wmax in pltwall
          if(bad_args(stack_req)) return
          values(15)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(15)
        case('t_pt')    ! total flight time of particles in caltraj
          if(bad_args(stack_req)) return
          values(16)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(16)
        case('rshi')    ! relative R shift to be applied to chords
          if(bad_args(stack_req)) return
          values(17)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(17)
          rshift=values(17)
        case('zshi')    ! relative Z shift to be applied to chords
          if(bad_args(stack_req)) return
          values(18)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(18)
          zshift=values(18)
        case('r0')      ! used to specify the R position of an ellipse in chord.F
          if(bad_args(stack_req)) return
          values(19)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(19)
          r0=values(19)
        case('z0')      ! used to specify the Z position of an ellipse in chord.F
          if(bad_args(stack_req)) return
          values(20)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(20)
          z0=values(20)
        case('a0')      ! used to specify the R radius of an ellipse in chord.F
          if(bad_args(stack_req)) return
          values(21)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(21)
          a0=values(21)
        case('b0')      ! used to specify the Z radius of an ellipse in chord.F
          if(bad_args(stack_req)) return
          values(22)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(22)
          b0=values(22)
        case('tt0')     ! used to specify the value within the ellipse in chord.F
          if(bad_args(stack_req)) return
          values(23)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(23)
          tt0=values(23)
        case('trgxmin')    ! Minimum x for target plotting in the computational domain
          if(bad_args(stack_req)) return
          values(24)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(24)
        case('trgxmax')    ! Maximum x for target plotting in the computational domain
          if(bad_args(stack_req)) return
          values(25)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(25)
        case('trgymin')    ! Minimum y for target plotting in the computational domain
          if(bad_args(stack_req)) return
          values(26)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(26)
        case('trgymax')    ! Maximum y for target plotting in the computational domain
          if(bad_args(stack_req)) return
          values(27)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(27)
        case('trgrmin')    ! Minimum r for target plotting in the physical domain
          if(bad_args(stack_req)) return
          values(28)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(28)
        case('trgrmax')    ! Maximum r for target plotting in the physical domain
          if(bad_args(stack_req)) return
          values(29)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(29)
        case('trgzmin')    ! Minimum z for target plotting in the physical domain
          if(bad_args(stack_req)) return
          values(30)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(30)
        case('trgzmax')    ! Maximum z for target plotting in the physical domain
          if(bad_args(stack_req)) return
          values(31)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(31)
        case('trgfmin')    ! Minimum function value for target plotting
          if(bad_args(stack_req)) return
          values(32)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(32)
        case('trgfmax')    ! Maximum function value for target plotting
          if(bad_args(stack_req)) return
          values(33)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(33)
        case('trgomin')    ! Minimum ordinate for 2-d target plots
          if(bad_args(stack_req)) return
          values(34)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(34)
        case('trgomax')    ! Maximum ordinate for 2-d target plots
          if(bad_args(stack_req)) return
          values(35)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(35)
        case('smin')       ! Minimum ordinate for wall plots
          if(bad_args(stack_req)) return
          values(36)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(36)
        case('smax')       ! Maximum ordinate for wall plots
          if(bad_args(stack_req)) return
          values(37)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(37)
        case default NVALUE
          cflag=.false.
      END SELECT NVALUE

      return

      end subroutine select_parameters

      subroutine select_command

      cflag=.true.
      stack_req=(/0,0,1,0,0,0/)
      SVALUE: select case(cval)
        case('extl')    ! extra label
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),stmp)
          call clear_extralabel(stmp)
        case('cltb')    ! colour table
          if(bad_args(stack_req)) return
          call set_color_map(spop(sstack,isstack))
        case('chor')    ! chord file
          if(bad_args(stack_req)) return
          call process_chords(spop(sstack,isstack))
        case('chpl')    ! plot chords
          if(bad_args(stack_req)) return
          call plot_chords(spop(sstack,isstack))
        case('cmd')     ! command file
          if(bad_args(stack_req)) return
          lu_pos=lu_pos+1
          if(lu_pos .gt. 10) then
            write(*,*) 'More than 10 nested files  --- infinite loop?'
            stop
          endif
          save_textline(lu_pos)=textline
          save_ipos(lu_pos)=ipos
          save_lu(lu_pos)=lu
          if(lu .eq. 5) then
            lu=50
          else
            lu=lu+1
          endif
          if(.not.open_file(lu,spop(sstack,isstack),
     1           'data.local/commands|data/commands')) then
            write(*,*) 'cmd failed top open file'
            lu=save_lu(lu_pos)
            lu_pos=lu_pos-1
          else
            write(*,*) 'reading from lu ',lu
            call readline(textline,lu,ipos)
          end if
        case('cmdcdf')     ! read te, ti, ne, na, dab2, dmb2 from the b2movies cdf file and loop calling a command file for each frame
          if(bad_args(stack_req)) return
          lu_pos=lu_pos+1
          if(lu_pos .gt. 10) then
            write(*,*) 'More than 10 nested files  --- infinite loop?'
            stop
          endif
          save_textline(lu_pos)=textline
          save_ipos(lu_pos)=ipos
          save_lu(lu_pos)=lu
          if(lu .eq. 5) then
            lu=50
          else
            lu=lu+1
          endif
          lucdf=lu
          if(.not.open_file(lu,spop(sstack,isstack),
     1           'data.local/commands|data/commands')) then
            write(*,*) 'cmd failed top open file'
            lu=save_lu(lu_pos)
            lu_pos=lu_pos-1
          else
#ifndef NO_CDF
            iret=nf_open('../b2movies.nc',NF_NOWRITE,ncidcdf)
            call check_cdf_status(iret)
            iret=nf_inq_varid(ncidcdf,'ntstep',varid)
            call check_cdf_status(iret)
            index(1)=1
            iret=nf_get_var1_int(ncidcdf,varid,index,ntstepcdf)
            call check_cdf_status(iret)
            write(*,*) ntstepcdf,' slices'
            istepcdf=1
            if(istepcdf.le.ntstepcdf) then
              call get_from_movie_cdf(ncidcdf,istepcdf)
              write(*,*) 'reading from lu ',lu,' for step ',istepcdf
              call readline(textline,lu,ipos)
            else
              lu=save_lu(lu_pos)
              lu_pos=lu_pos-1
            endif
#else
            write(*,*) 'Not compiled with CDF'
            lu=save_lu(lu_pos)
            lu_pos=lu_pos-1
#endif
          end if
        case('chbx')    ! boxcar average of chords
          if(bad_args(stack_req)) return
          call process_chords_boxcar(spop(sstack,isstack))
        case('chvl')    ! value at positions 
          if(bad_args(stack_req)) return
          call process_chord_slice(spop(sstack,isstack))
        case('flux')    ! call WOSlib flux routine
          if(bad_args(stack_req)) return
          call flux(spop(sstack,isstack),nx,ny,
     1           values(5),values(6),values(7),values(8))
        case('chmp')    ! Plots a contour plot of the chord contribution to the mesh element
          if(bad_args(stack_req)) return
          call map_chords(spop(sstack,isstack))
        case('cam')     ! camera plot
          if(bad_args(stack_req)) return
          call process_camera(spop(sstack,isstack))
        case('movi')    ! movie
          if(bad_args(stack_req)) return
          moviename=spop(sstack,isstack)
        case('sonn')    ! something to do with Sonnet (no longer working)
          if(bad_args(stack_req)) return
c          sonnetfile=spop(sstack,isstack)
c          write(*,*) 'Reading a new geometry'
c          iuserinit=0
c          read_xdr_geometry=.true.
c          call b2igeo
c          write(*,*) 'Read a new geometry'
          stmp=spop(sstack,isstack)
          write(*,*) 'Read a new Sonnet geometry : No longer supported!'
        case('saw')     ! reads external data onto the matrix stack
          stack_req=(/2,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          arg=spop(sstack,isstack)
          ny=ipop(istack,iistack)-2
          nx=ipop(istack,iistack)-2
          nz=1
          open(99,file=arg,status='old')
          do ix=-1,nx
            do iy=-1,ny
              read(99,*) cr(ix,iy),cz(ix,iy),
     1     crx(ix,iy,1),cry(ix,iy,1),
     2     crx(ix,iy,3),cry(ix,iy,3),
     3     crx(ix,iy,0),cry(ix,iy,0),
     4     crx(ix,iy,2),cry(ix,iy,2),
     5     mach(ix,iy,0)
            enddo
          enddo
          do IX = -1,NX
            XX(IX,2) = IX
            XX(IX,1) = IX-0.5_R8
          end do
          do IY = -1, NY
            YY(IY,2) = IY
            YY(IY,1) = IY-0.5_R8
          end do
          values(1)=minval(xx)                             	 !xmin
          values(2)=maxval(xx)+0.5_R8                      	 !xmax
          values(3)=minval(yy)                             	 !ymin
          values(4)=maxval(yy)+0.5_R8                      	 !ymax
          values(5)=minval(crx)                            	 !rmin
          values(6)=maxval(crx)                            	 !rmax
          values(7)=minval(cry)                            	 !zmin
          values(8)=maxval(cry)                            	 !zmax
          write(*,*) 'values',(i,values(i),i=1,8)
          call copy(mach,1.0_R8,nz,textbl)
          graphlabel(iwrk)=arg
        case('oson')    ! output geometry in Sonnet format
          if(bad_args(stack_req)) return
          arg=spop(sstack,isstack)
          open(99,file=arg)
          icell=0
          write(99,*)
          write(99,'(3x,a/)') 'Element output:'
          if(values(19).eq.0.0_R8 .or. values(22).eq.0.0_R8) then
            write(*,*) 'Using default values for B0 and R0 !'
            b0r0=abs(cr(-1,-1)*bb(-1,-1,2))
          else
            b0r0=values(22)*values(19)
          endif
          write(99,*) '   R*Btor =',b0r0
          write(99,'(/3x,2a)')
     1           '======================================',
     2           '=================================================='
          do iy=-1,ny
            do ix=-1,nx
              write(99,10001) icell,ix+1,iy+1,
     1         crx(ix,iy,2),cry(ix,iy,2),crx(ix,iy,3),cry(ix,iy,3)
10001         format(3x,'Element',i5,' = (',i4,',',i4,'): (',
     .                 1pe17.10,',',1pe17.10,')',
     .                 6x,'(',1pe17.10,',',1pe17.10,')')
              write(99,10002)
     1         sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3),
     2         cr(ix,iy),cz(ix,iy)
10002         format(3x,'Field ratio  = ',1pe17.10,13x,
     .                 '(',1pe17.10,',',1pe17.10,')')
              write(99,10003)
     1         crx(ix,iy,0),cry(ix,iy,0),crx(ix,iy,1),cry(ix,iy,1)
10003         format(
     .                 t30,'(',1pe17.10,',',1pe17.10,')',
     .                 6x,'(',1pe17.10,',',1pe17.10,')')
              write(99,10005)
10005         format(3x,
     1         '--------------------------------------------------',
     2         '--------------------------------------')
              icell=icell+1
            enddo
          enddo
          close(99)
        case('uedgegrd')    ! output geometry in UEDGE format (string + real arg)
          stack_req=(/0,1,1,0,0,0/)
          if(bad_args(stack_req)) return
          arg=spop(sstack,isstack)
          t1=rpop(rstack,irstack)
          open(99,file=arg)
c-dpc still in testing phase
          write(99,'(5i4)') nx,ny,leftcut(1),rightcut(1),topcut(1)
          write(99,'()')       	 ! rm
          write(99,'(1p,3e23.15)') ((cr(ix,iy),ix=-1,nx),iy=-1,ny),
     1     ((crx(ix,iy,0),ix=-1,nx),iy=-1,ny),
     1     ((crx(ix,iy,1),ix=-1,nx),iy=-1,ny),
     1     ((crx(ix,iy,2),ix=-1,nx),iy=-1,ny),
     1     ((crx(ix,iy,3),ix=-1,nx),iy=-1,ny)
          write(99,'()')       	 ! zm
          write(99,'(1p,3e23.15)') ((cz(ix,iy)+t1,ix=-1,nx),iy=-1,ny),
     1     ((cry(ix,iy,0)+t1,ix=-1,nx),iy=-1,ny),
     1     ((cry(ix,iy,1)+t1,ix=-1,nx),iy=-1,ny),
     1     ((cry(ix,iy,2)+t1,ix=-1,nx),iy=-1,ny),
     1     ((cry(ix,iy,3)+t1,ix=-1,nx),iy=-1,ny)
          write(99,'()')       	 ! psi
          uedge_c=0.0_R8
          uedge_v=0.0_R8
          write(99,'(1p,3e23.15)') ((uedge_c(ix,iy),ix=-1,nx),iy=-1,ny),
     1     (((uedge_v(ix,iy,in),ix=-1,nx),iy=-1,ny),in=1,4)
          write(99,'()')       	 ! br
          uedge_c=0.0_R8
          uedge_v=0.0_R8
          write(99,'(1p,3e23.15)') ((uedge_c(ix,iy),ix=-1,nx),iy=-1,ny),
     1     (((uedge_v(ix,iy,in),ix=-1,nx),iy=-1,ny),in=1,4)
          write(99,'()')       	 ! bz
          uedge_c=0.0_R8
          uedge_v=0.0_R8
          write(99,'(1p,3e23.15)') ((uedge_c(ix,iy),ix=-1,nx),iy=-1,ny),
     1     (((uedge_v(ix,iy,in),ix=-1,nx),iy=-1,ny),in=1,4)
          write(99,'()')       	 ! bpol
          uedge_c(:,:)=bb(:,:,0)
          write(*,*) 'bpol_c: min & max ',
     1     minval(uedge_c),maxval(uedge_c)
          uedge_v=0.0_R8
          call my_intcor(nx,ny,vol,uedge_c,corners)
          write(*,*) 'bpol_v: min & max ',
     1     minval(corners),maxval(corners)
          do iy=-1,ny
            do ix=-1,nx
              if (b2plot_debug.ge.12) write(*,*) 'ix,iy = ',ix,iy
              uedge_v(ix,iy,1)=corners(ix,iy)    ! bottom left
              ixr=rightix(ix,iy)
              iyt=topiy(ix,iy)
              if (b2plot_debug.ge.12) write(*,*) 'ixr,iy = ',ixr,iy
              uedge_v(ix,iy,2)=corners(ixr,iy)   ! bottom right
              if (b2plot_debug.ge.12) write(*,*) 'ix,iyt = ',ix,iyt
              uedge_v(ix,iy,3)=corners(ix,iyt)   ! top left
              if(iyt.le.ny) then
                ixr=rightix(ix,iyt)
              else if(ixr.le.nx) then
                iyt=topiy(ixr,iy)
              else
                ixr=ix+1
                iyt=iy+1
              endif
              if (b2plot_debug.ge.12) write(*,*) 'ixr,iyt = ',ixr,iyt
              uedge_v(ix,iy,4)=corners(ixr,iyt)  ! top right
            enddo
          enddo
          write(*,*) 'bpol_v: min & max ',
     1     minval(uedge_v),maxval(uedge_v)
          write(99,'(1p,3e23.15)') ((uedge_c(ix,iy),ix=-1,nx),iy=-1,ny),
     1     (((uedge_v(ix,iy,in),ix=-1,nx),iy=-1,ny),in=1,4)
          write(99,'()')       	 ! bphi
          uedge_c(:,:)=bb(:,:,2)
          write(*,*) 'bphi_c: min & max ',
     1     minval(uedge_c),maxval(uedge_c)
          uedge_v=0.0_R8
          call my_intcor(nx,ny,vol,uedge_c,corners)
          write(*,*) 'bphi_v: min & max ',
     1     minval(corners),maxval(corners)
          do iy=-1,ny
            do ix=-1,nx
              if (b2plot_debug.ge.12) write(*,*) 'ix,iy = ',ix,iy
              uedge_v(ix,iy,1)=corners(ix,iy)    ! bottom left
              ixr=rightix(ix,iy)
              iyt=topiy(ix,iy)
              if (b2plot_debug.ge.12) write(*,*) 'ixr,iy = ',ixr,iy
              uedge_v(ix,iy,2)=corners(ixr,iy)   ! bottom right
              if (b2plot_debug.ge.12) write(*,*) 'ix,iyt = ',ix,iyt
              uedge_v(ix,iy,3)=corners(ix,iyt)   ! top left
              if(iyt.le.ny) then
                ixr=rightix(ix,iyt)
              else if(ixr.le.nx) then
                iyt=topiy(ixr,iy)
              else
                ixr=ix+1
                iyt=iy+1
              endif
              if (b2plot_debug.ge.12) write(*,*) 'ixr,iyt = ',ixr,iyt
              uedge_v(ix,iy,4)=corners(ixr,iyt)  ! top right
            enddo
          enddo
          write(*,*) 'bphi_v: min & max ',
     1     minval(uedge_v),maxval(uedge_v)
          write(99,'(1p,3e23.15)') ((uedge_c(ix,iy),ix=-1,nx),iy=-1,ny),
     1     (((uedge_v(ix,iy,in),ix=-1,nx),iy=-1,ny),in=1,4)
          write(99,'()')       	 ! b
          uedge_c(:,:)=bb(:,:,3)
          write(*,*) 'b_c: min & max ',
     1     minval(uedge_c),maxval(uedge_c)
          uedge_v=0.0_R8
          call my_intcor(nx,ny,vol,uedge_c,corners)
          write(*,*) 'b_v: min & max ',
     1     minval(corners),maxval(corners)
          do iy=-1,ny
            do ix=-1,nx
              if (b2plot_debug.ge.12) write(*,*) 'ix,iy = ',ix,iy
              uedge_v(ix,iy,1)=corners(ix,iy)    ! bottom left
              ixr=rightix(ix,iy)
              iyt=topiy(ix,iy)
              if (b2plot_debug.ge.12) write(*,*) 'ixr,iy = ',ixr,iy
              uedge_v(ix,iy,2)=corners(ixr,iy)   ! bottom right
              if (b2plot_debug.ge.12) write(*,*) 'ix,iyt = ',ix,iyt
              uedge_v(ix,iy,3)=corners(ix,iyt)   ! top left
              if(iyt.le.ny) then
                ixr=rightix(ix,iyt)
              else if(ixr.le.nx) then
                iyt=topiy(ixr,iy)
              else
                ixr=ix+1
                iyt=iy+1
              endif
              if (b2plot_debug.ge.12) write(*,*) 'ixr,iyt = ',ixr,iyt
              uedge_v(ix,iy,4)=corners(ixr,iyt)  ! top right
            enddo
          enddo
          write(*,*) 'b_v: min & max ',
     1     minval(uedge_v),maxval(uedge_v)
          write(99,'(1p,3e23.15)') ((uedge_c(ix,iy),ix=-1,nx),iy=-1,ny),
     1     (((uedge_v(ix,iy,in),ix=-1,nx),iy=-1,ny),in=1,4)
          write(99,'(a60)') 'output from b2pl.exe'
          close(99)
        case('expf')    ! add experimental data to plot
          if(bad_args(stack_req)) return
          expfile=spop(sstack,isstack)
          call plot_expfile
        case('outp')    ! prints the value of a specific element
          stack_req=(/2,0,1,1,0,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          iy=ipop(istack,iistack)
          ix=ipop(istack,iistack)
          do is=0,ndim(3,iwrk)-1
            if(ndim(3,iwrk) .eq. 1) then
              write(*,*) stmp(1:lnblnk(stmp)),
     1                 work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              if(writing_to_file) then
                write(nwtf,*) stmp(1:lnblnk(stmp)),                   !srv 21.09.11
     1                    work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              endif
            else
              write(*,*) stmp(1:lnblnk(stmp)),is,
     1                 work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              if(writing_to_file) then
                write(nwtf,*) stmp(1:lnblnk(stmp)),is,                 !srv 21.09.11
     1                    work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              endif
            endif
          enddo
        case('pdev')    ! PLPLOT device name
          if(bad_args(stack_req)) return
          call set_device(spop(sstack,isstack))
        case('mprt')    ! output the top of matrix stack to a file
          stack_req=(/0,0,1,1,0,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          open(99,file=stmp)
          do is=0,ndim(3,iwrk)-1
            write(99,'(''#'',1x,a,1x,a)')
     1       graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))),
     2       textwk(is,iwrk)(1:lnblnk(textwk(is,iwrk)))
            do iy=-1,ny
              do ix=-1,nx
                write(99,*) ix,iy,is,work(ix,iy,is,iwrk)
              enddo
            enddo
          enddo
          close(99)
        case('mtab')    ! output a table from the top of matrix stack to a file
          stack_req=(/4,0,1,1,0,0/)
          if(bad_args(stack_req)) return
          iy2=ipop(istack,iistack)
          iy1=ipop(istack,iistack)
          if (iy1.eq.0 .and. iy2.eq.0) then
            iy1 = -1
            iy2 = ny
          end if
          ix2=ipop(istack,iistack)
          ix1=ipop(istack,iistack)
          if (ix1.eq.0 .and. ix2.eq.0) then
            ix1 = -1
            ix2 = nx
          end if
          if (ix1 .lt. -1 .or. ix2 .gt. nx .or.
     &        iy1 .lt. -1 .or. iy2 .gt. ny .or.
     1        ix1 .gt. ix2 .or. iy1 .gt. iy2) then
            write(*,*) 'Parameters supplied to mtab were invalid - ',
     1                 'valid parameters are:'
            write(*,*) ' -1 <= ',ix1,' <= ',ix2,' <= ',nx
            write(*,*) ' -1 <= ',iy1,' <= ',iy2,' <= ',ny
            return
          end if
          stmp=spop(sstack,isstack)
          nz=ndim(3,iwrk)
          open(99,file=stmp)
          write(hlp_frm,'(a,i2,a)') '(a1,1x,a,1x,',nz,'a13)'
          write(99,hlp_frm) '#',
     1       graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))),
     2       (textwk(is,iwrk)(1:lnblnk(textwk(is,iwrk))),is=0,nz-1)
          write(hlp_frm,'(a,i2,a)') '(2i4,1p,',nz,'(1e12.5,1x))'
          do iy=iy1,iy2
            do ix=ix1,ix2
              write(99,hlp_frm) ix,iy,work(ix,iy,0:nz-1,iwrk)
            enddo
          enddo
          close(99)
        case('amju','amjuel')    ! accesses Eirene type atomic physics data
          stack_req=(/0,0,4,0,0,0/)
          if(bad_args(stack_req)) return
          crc=spop(sstack,isstack)
          reac=spop(sstack,isstack)
          h123=spop(sstack,isstack)
          filnam=spop(sstack,isstack)
!#ifdef USE_EIRENE_NOT_YET_READY
          call get_amjuel_data(filnam,h123,reac,crc,
     1     nxd,nyd,nx,ny,ne,te,mach)
!#endif
          call copy(mach,1.0_R8,1,textbl)
          graphlabel(iwrk)=filnam//h123//reac//crc
        case('scat')    ! concatenate two strings on top of string stack
          stack_req=(/0,0,2,0,0,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          stmp=trim(spop(sstack,isstack)) // stmp
          call spush(sstack,isstack,stmp)
        case('sdrop')   ! drop top of string stack
          stmp=spop(sstack,isstack)
        case('sdup')    ! duplicate top of string stack
          stmp=sstack(isstack)
          call spush(sstack,isstack,stmp)
        case('sswap')   ! swap top two elements of string stack
          stack_req=(/0,0,2,0,0,0/)
          if(bad_args(stack_req)) return
          stmp1=spop(sstack,isstack)
          stmp2=spop(sstack,isstack)
          call spush(sstack,isstack,stmp1)
          call spush(sstack,isstack,stmp2)
        case('sover')   ! str2 str1 --> str2 str1 str2
          stack_req=(/0,0,2,0,0,0/)
          if(bad_args(stack_req)) return
          stmp=sstack(isstack-1)
          call spush(sstack,isstack,stmp)
        case('sprt')    ! output the string stack
          stack_req=(/0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          do is = isstack, 1, -1
            write(*,'(a)') trim(sstack(is))
            if(writing_to_file) then
              write(nwtf,'(a)') trim(sstack(is))                   !srv 21.09.11
            endif
          enddo
        case('omovie')      ! old movie option
          stack_req=(/2,0,2,0,0,0/)
          if(bad_args(stack_req)) return
          varname=spop(sstack,isstack)
          filename=spop(sstack,isstack)
          iend=ipop(istack,iistack)
          istart=ipop(istack,iistack)
          itmp=0
#ifndef NO_CDF
          call get_cdf_ntstep(filename,itmp)
          write(*,*) 'MOVIE: ',itmp
          if(iend.lt.1) iend=itmp
          if(iwrk .eq. wklng) call dropstackbottom()
          iwrk=iwrk+1
          ndim(1,iwrk)=nx
          ndim(2,iwrk)=ny
          ndim(3,iwrk)=1
          textwk(:,iwrk)=' '
          graphlabel(iwrk)=' '
          do is=max(istart,1),min(iend,itmp)
            call get_cdf(filename, is, movietime(1),
     1       varname, mach, nx, ny)
            write(*,*) 'MOVIE: ',movietime(1)
            write(stmp,'(f13.9,2x,i5,''/'',i5)') movietime(1), is, itmp
            description=trim(filename) // ' ' // trim(varname) //
     1       ' ' // stmp(1:26)
            write(*,*) 'MOVIE: ',trim(description)
            work(:,:,0:0,iwrk)=mach(:,:,0:0)
            graphlabel(iwrk)=trim(description)
            call surf(min(iend,itmp)-max(istart,1)+1)
          enddo
          iwrk=iwrk-1
#else
          write(*,*) 'Compile b2plot without NO_CDF option for movies!'
#endif
        case('movie','nmovie')    ! new movie option
          stack_req=(/2,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          varname=spop(sstack,isstack)
          filename='../b2movies.nc'
          iend=ipop(istack,iistack)
          istart=ipop(istack,iistack)
          itmp=0
#ifndef NO_CDF
          call get_cdf_ntstep(filename,itmp)
          write(*,*) 'MOVIE: ',itmp
          if(iend.lt.1) iend=itmp
          if(iwrk .eq. wklng) call dropstackbottom()
          iwrk=iwrk+1
          ndim(1,iwrk)=nx
          ndim(2,iwrk)=ny
          ndim(3,iwrk)=1
          textwk(:,iwrk)=' '
          graphlabel(iwrk)=' '
          do is=max(istart,1),min(iend,itmp)
            call get_cdf(filename, is, movietime(1),
     1       varname, mach, nx, ny)
            write(*,*) 'MOVIE: ',movietime(1)
            write(stmp,'(f13.9,2x,i5,''/'',i5)') movietime(1), is, itmp
            description=trim(filename) // ' ' // trim(varname) //
     1       ' ' // stmp(1:26)
            write(*,*) 'MOVIE: ',trim(description)
            work(:,:,0:0,iwrk)=mach(:,:,0:0)
            graphlabel(iwrk)=trim(description)
            call surf(min(iend,itmp)-max(istart,1)+1)
          enddo
          iwrk=iwrk-1
#else
          write(*,*) 'Compile b2plot without NO_CDF option for movies!'
#endif
        case('wallblbl')  ! set wall bar label
          if(bad_args(stack_req)) return
          call set_wallbarlabel(spop(sstack,isstack))
        case('trgblbl')  ! set target bar label
          if(bad_args(stack_req)) return
          call set_targetbarlabel(spop(sstack,isstack))
        case('glbl')  ! set graph label
          stack_req=(/0,0,1,1,0,0/)
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),graphlabel(iwrk))
        case('wlbl')  ! set wall graph label
          stack_req=(/0,0,1,0,1,0/)
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),wall_label(iwwrk))
        case('tlbl')  ! set target graph label
          stack_req=(/0,0,1,0,0,1/)
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),target_label(itwrk))
        case('slbl')  ! set species label
          stack_req=(/1,0,1,1,0,0/)
          if(bad_args(stack_req)) return
          is1=ipop(istack,iistack)
          if (0.gt.is1.or.is1.ge.ndim(3,iwrk)) then
            write(*,*) 'Faulty parameter is = ',is1,' for slbl!'
            return
          endif
          textwk(is1,iwrk)=spop(sstack,isstack)
        case('wslbl')  ! set wall species label
          stack_req=(/1,0,1,0,1,0/)
          if(bad_args(stack_req)) return
          is1=ipop(istack,iistack)
          if (0.gt.is1.or.is1.ge.nwdim(2,iwwrk)) then
            write(*,*) 'Faulty parameter is = ',is1,' for wslbl!'
            return
          endif
          textwwk(is1,iwwrk)=spop(sstack,isstack)
        case('tranprof')  ! output transport profile
          stack_req=(/0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          open(99,file='b2.transport.inputfile')
          write(99,*) '&TRANSPORT'
          r_ref=(crx(jxa,jsep,2)+crx(jxa,jsep,3))/2.0_R8
          z_ref=(cry(jxa,jsep,2)+cry(jxa,jsep,3))/2.0_R8
          do iy=-1,jsep
            s_ref(iy)=-sqrt((cr(jxa,iy)-r_ref)**2+
     1       (cz(jxa,iy)-z_ref)**2)
          enddo
          do iy=jsep+1,ny
            s_ref(iy)=sqrt((cr(jxa,iy)-r_ref)**2+
     1       (cz(jxa,iy)-z_ref)**2)
          enddo
          s_ref(-1)=1.1_R8 * s_ref(-1)
          s_ref(ny)=1.1_R8 * s_ref(ny)
          do is=0,ns-1
            if(.not.is_neutral(is)) then
              do iy=-1,ny
                profile(iy)=dna0(jxa,iy,is)
              enddo
              call write_crazy(99,ny+2,1,is,s_ref,profile)
            endif
          enddo
          do is=0,ns-1
            if(.not.is_neutral(is)) then
              do iy=-1,ny
                profile(iy)=hcib(jxa,iy,is)/na(jxa,iy,is)
              enddo
              call write_crazy(99,ny+2,3,is,s_ref,profile)
            endif
          enddo
          do iy=-1,ny
            profile(iy)=hce0(jxa,iy)/ne(jxa,iy)
          enddo
          call write_crazy(99,ny+2,4,0,s_ref,profile)
          write(99,*) '/'
          close(99)
        case default SVALUE
          cflag=.false.
      END SELECT SVALUE

      if(cflag.and.b2plot_debug.ge.1) write(*,*) trim(cval),' : Done.'
      return

      end subroutine select_command

      subroutine select_wall_stack

      cflag=.true.
      WVALUE: select case(cval)
        case('wallarea')
          description='Wall area (m^{2})'
          nz=1
          call wall_copy(wall_area,1.0_R8,nz,textbl)
        case('ei')
          description='Incident energy (eV)'
          nz=nsw
          call wall_copy(incident_energy,1.0_R8,nz,textwl)
        case('wallflux')
          description='Wall particle flux (s-1)'
          nz=nsw
          call wall_copy(wall_flux,1.0_R8,nz,textpl)
        case('wallhflx')
          description='Wall heat flux (W/m^{2})'
          nz=1
          call wall_copy(wall_heatflux,1.0_R8,nz,textbl)
        case('r1')
          description='Wall start R position (m)'
          nz=1
          call wall_copy(r1,1.0_R8,nz,textbl)
        case('r2')
          description='Wall end R position (m)'
          nz=1
          call wall_copy(r2,1.0_R8,nz,textbl)
        case('z1')
          description='Wall start Z position (m)'
          nz=1
          call wall_copy(z1,1.0_R8,nz,textbl)
        case('z2')
          description='Wall end Z position (m)'
          nz=1
          call wall_copy(z2,1.0_R8,nz,textbl)
        case('walltemp')
          description='Wall surface temperature (K)'
          nz=1
          call wall_copy(wall_temp,1.0_R8,nz,textbl)
        case('plathick')
          description='Plate thickness (m)'
          nz=1
          call wall_copy(plate_thickness,1.0_R8,nz,textbl)
        case('coathick')
          description='Coating thickness (m)'
          nz=1
          call wall_copy(coating_thickness,1.0_R8,nz,textbl)
        case('cooltemp')
          description='Backplate temperature (K)'
          nz=1
          call wall_copy(backplate_temp,1.0_R8,nz,textbl)
        case('vaprtemp')
          description='Vaporisation temperature (K)'
          nz=1
          call wall_copy(matsurf_vapor_temp,1.0_R8,nz,textbl)
        case('north')
          description='North boundary'
          nz=1
          do ibnd = 1, nbnd
            if(topiy(imapx(ibnd),imapy(ibnd)).eq.ny+1) then
              data(ibnd) = 1.0_R8
            else
              data(ibnd) = 0.0_R8
            endif
          enddo
          call wall_copy(data,1.0_R8,nz,textbl)
        case('south')
          description='South boundary'
          nz=1
          do ibnd = 1, nbnd
            if(bottomiy(imapx(ibnd),imapy(ibnd)).eq.-2) then
              data(ibnd) = 1.0_R8
            else
              data(ibnd) = 0.0_R8
            endif
          enddo
          call wall_copy(data,1.0_R8,nz,textbl)
        case('west')
          description='West boundary'
          nz=1
          do ibnd = 1, nbnd
            if(leftix(imapx(ibnd),imapy(ibnd)).eq.-2) then
              data(ibnd) = 1.0_R8
            else
              data(ibnd) = 0.0_R8
            endif
          enddo
          call wall_copy(data,1.0_R8,nz,textbl)
        case('east')
          description='East boundary'
          nz=1
          do ibnd = 1, nbnd
            if(rightix(imapx(ibnd),imapy(ibnd)).eq.nx+1) then
              data(ibnd) = 1.0_R8
            else
              data(ibnd) = 0.0_R8
            endif
          enddo
          call wall_copy(data,1.0_R8,nz,textbl)
        case('nseg')
          description='Wall segment'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(abs(istack(iistack)).gt.nsegments) then
            write(*,*) 'No such wall segment'
            return
          endif
          if(istack(iistack) .eq. 0) then
            do i = 1, nsegments
              data(wall_begin(i):wall_end(i)) = real(i)
            enddo
          elseif(istack(iistack) .gt. 0) then
            data(wall_begin(1):wall_end(istack(iistack)-1)) = 0.0_R8
            data(wall_begin(istack(iistack)):
     .           wall_end(istack(iistack))) = 1.0_R8
            data(wall_begin(istack(iistack)+1):
     .           wall_end(nsegments)) = 0.0_R8
          else
            data(wall_begin(1):
     .           wall_end(abs(istack(iistack))-1)) = 1.0_R8
            data(wall_begin(abs(istack(iistack))):
     .           wall_end(abs(istack(iistack)))) = 0.0_R8
            data(wall_begin(abs(istack(iistack))+1):
     .           wall_end(nsegments)) = 1.0_R8
          endif
          write(*,*) '"nseg" option with ',ipop(istack,iistack)
          call wall_copy(data,1.0_R8,nz,textbl)
        case ('iwlt')
          description='iwlt'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          data(1:istack(iistack)-1) = 1.0_R8
          data(istack(iistack):nbnd) = 0.0_R8
          write(*,*) '"iwlt" option with ',ipop(istack,iistack)
          call wall_copy(data,1.0_R8,nz,textbl)
        case ('iweq')
          description='iweq'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          data(1:nbnd) = 0.0_R8
          data(istack(iistack)) = 1.0_R8
          write(*,*) '"iweq" option with ',ipop(istack,iistack)
          call wall_copy(data,1.0_R8,nz,textbl)
        case ('iwgt')
          description='iwgt'
          nz=1
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          data(1:istack(iistack)) = 0.0_R8
          data(istack(iistack)+1:nbnd) = 1.0_R8
          write(*,*) '"iwgt" option with ',ipop(istack,iistack)
          call wall_copy(data,1.0_R8,nz,textbl)
        case ('iswlt')
          description='iswlt'
          nz=nsw
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wall_tmp(1:nbnd,0:istack(iistack)-1) = 1.0_R8
          wall_tmp(1:nbnd,istack(iistack):ns-1) = 0.0_R8
          write(*,*) '"iswlt" option with ',ipop(istack,iistack)
          call wall_copy(wall_tmp,1.0_R8,nz,textwl)
        case ('isweq')
          description='isweq'
          nz=nsw
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wall_tmp(1:nbnd,0:ns-1) = 0.0_R8
          wall_tmp(1:nbnd,istack(iistack)) = 1.0_R8
          write(*,*) '"isweq" option with ',ipop(istack,iistack)
          call wall_copy(wall_tmp,1.0_R8,nz,textwl)
        case ('iswgt')
          description='iswgt'
          nz=nsw
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wall_tmp(1:nbnd,0:istack(iistack)) = 0.0_R8
          wall_tmp(1:nbnd,istack(iistack)+1:ns-1) = 1.0_R8
          write(*,*) '"iswgt" option with ',ipop(istack,iistack)
          call wall_copy(wall_tmp,1.0_R8,nz,textwl)
        case('ibnd')
          description='ibnd'
          nz=1
          do i = 1, nbnd
            data(i) = real(i)
          enddo
          call wall_copy(data,1.0_R8,nz,textbl)
        case('walldist')
          description='Distance along wall (m)'
          nz=1
          do i = 1, nbnd
            if (i.eq.1) then
              data(i) = 0.5_R8*
     .         sqrt((r2(i)-r1(i))**2+(z2(i)-z1(i))**2)
            else
              data(i) = data(i-1) + 0.5_R8*(
     .         sqrt((r2(i-1)-r1(i-1))**2+(z2(i-1)-z1(i-1))**2)+
     .         sqrt((r2(i)-r1(i))**2+(z2(i)-z1(i))**2))+
     .         sqrt((r1(i)-r2(i-1))**2+(z1(i)-z2(i-1))**2)
            endif
          enddo
          call wall_copy(data,1.0_R8,nz,textbl)
       case('wallngth')
          description='Wall length (m)'
          nz=1
          do i = 1, nbnd
            data(i) = sqrt((r2(i)-r1(i))**2+(z2(i)-z1(i))**2)
          enddo
          call wall_copy(data,1.0_R8,nz,textbl)
        case('sput')
          if (sput_dst.lt.0 .or. sput_dst.ge.ns) return
          description='Physical sputtering yield of'//
     .     trim(textpl(sput_dst))//' by'
          nz=nsw
          call wall_copy(physical_sputtering,1.0_R8,nz,textwl)
        case('sput2')
          if (sput_dst.lt.0 .or. sput_dst.ge.ns) return
          description='Physical sputtering energy yield of'//
     .     trim(textpl(sput_dst))//' by'
          nz=nsw
          call wall_copy(physical_sputtering_energy,1.0_R8,nz,textwl)
        case('sputchem')
          description='Chemical sputtering yield'
          nz=nsw
          call wall_copy(chemical_sputtering,1.0_R8,nz,textwl)
        case('thermev')
          description='Thermal Evaporation Fluxes (/m2/s)'
          nz=nsw
          call wall_copy(thermal_evaporation,1.0_R8,nz,textwl)
        case('sputres')
          description='RES sputtering yield'
          nz=nsw
          call wall_copy(RES_sputtering,1.0_R8,nz,textwl)
        case('backscat')
          description='Backscattering fraction'
          nz=nsw
          call wall_copy(backscattering,1.0_R8,nz,textwl)
        case('backscae')
          description='Backscattering energy fraction'
          nz=nsw
          call wall_copy(backscattering_energy,1.0_R8,nz,textwl)
        case('erosion')
          description='Gross erosion'
          nz=ntargsp
          call wall_copy(erosion,1.0_R8,nz,target_code)
        case('deposit')
          description='Gross deposition'
          nz=ntargsp
          call wall_copy(deposition,1.0_R8,nz,target_code)
        case('ero_ml')
          description='Monolayer erosion'
          nz=ntargsp
          call wall_copy(monolayer_erosion,1.0_R8,nz,target_code)
        case('depo_ml')
          description='Monolayer deposition'
          nz=ntargsp
          call wall_copy(monolayer_deposition,1.0_R8,nz,target_code)
        case('time')
          description='Plate time (s)'
          nz=1
          call wall_copy(plate_time,1.0_R8,nz,textbl)
        case('imapx')
          description='imapx'
          nz=1
          data(1:nbnd)=real(imapx(1:nbnd))
          call wall_copy(data,1.0_R8,nz,textbl)
        case('imapy')
          description='imapy'
          nz=1
          data(1:nbnd)=real(imapy(1:nbnd))
          call wall_copy(data,1.0_R8,nz,textbl)
        case('kappa_a')
          description='Bulk Kappa A factor (1/AT+B) (m/W)'
          nz=1
          call wall_copy(matbulk_Kappa_KA,1.0_R8,nz,textbl)
        case('kappa_b')
          description='Bulk Kappa B factor (1/AT+B) (m.K/W)'
          nz=1
          call wall_copy(matbulk_Kappa_KB,1.0_R8,nz,textbl)
        case('kappa_a1')
          description='Bulk Kappa A1 polynomial factor (W/m/K)'
          nz=1
          call wall_copy(matbulk_Kappa_A1,1.0_R8,nz,textbl)
        case('kappa_a2')
          description='Bulk Kappa A2 polynomial factor (W/m/K^{2})'
          nz=1
          call wall_copy(matbulk_Kappa_A2,1.0_R8,nz,textbl)
        case('kappa_a3')
          description='Bulk Kappa A3 polynomial factor (W/m/K^{3})'
          nz=1
          call wall_copy(matbulk_Kappa_A3,1.0_R8,nz,textbl)
        case('kappa_a4')
          description='Bulk Kappa A4 polynomial factor (W/m/K^{4})'
          nz=1
          call wall_copy(matbulk_Kappa_A4,1.0_R8,nz,textbl)
        case('kappapa1')
          description='Bulk Kappa_perp A1 polynomial factor (W/m/K)'
          nz=1
          call wall_copy(matbulk_Kappa_perp_A1,1.0_R8,nz,textbl)
        case('kappapa2')
          description='Bulk Kappa_perp A2 polynomial factor (W/m/K^{2})'
          nz=1
          call wall_copy(matbulk_Kappa_perp_A2,1.0_R8,nz,textbl)
        case('kappapa3')
          description='Bulk Kappa_perp A3 polynomial factor (W/m/K^{3})'
          nz=1
          call wall_copy(matbulk_Kappa_perp_A3,1.0_R8,nz,textbl)
        case('kappapa4')
          description='Bulk Kappa_perp A4 polynomial factor (W/m/K^{4})'
          nz=1
          call wall_copy(matbulk_Kappa_perp_A4,1.0_R8,nz,textbl)
        case('tgdns')
          description='Bulk mass density (kg/m^{3})'
          nz=1
          call wall_copy(matbulk_target_density,1.0_R8,nz,textbl)
        case('cp_b1')
          description='Bulk heat capacity B1 factor (J/g/K)'
          nz=1
          call wall_copy(matbulk_Cp_B1,1.0_R8,nz,textbl)
        case('cp_b2')
          description='Bulk heat capacity B2 factor (J/g/K^{2})'
          nz=1
          call wall_copy(matbulk_Cp_B2,1.0_R8,nz,textbl)
        case('cp_b3')
          description='Bulk heat capacity B3 factor (J.K/g)'
          nz=1
          call wall_copy(matbulk_Cp_B3,1.0_R8,nz,textbl)
        case('ckappa_a')
          description='Coating Kappa A factor (1/AT+B) (m/W)'
          nz=1
          call wall_copy(coating_Kappa_KA,1.0_R8,nz,textbl)
        case('ckappa_b')
          description='Coating Kappa B factor (1/AT+B) (m.K/W)'
          nz=1
          call wall_copy(coating_Kappa_KB,1.0_R8,nz,textbl)
        case('ckappaa1')
          description='Coating Kappa A1 polynomial factor (W/m/K)'
          nz=1
          call wall_copy(coating_Kappa_A1,1.0_R8,nz,textbl)
        case('ckappaa2')
          description='Coating Kappa A2 polynomial factor (W/m/K^{2})'
          nz=1
          call wall_copy(coating_Kappa_A2,1.0_R8,nz,textbl)
        case('ckappaa3')
          description='Coating Kappa A3 polynomial factor (W/m/K^{3})'
          nz=1
          call wall_copy(coating_Kappa_A3,1.0_R8,nz,textbl)
        case('ckappaa4')
          description='Coating Kappa A4 polynomial factor (W/m/K^{4})'
          nz=1
          call wall_copy(coating_Kappa_A4,1.0_R8,nz,textbl)
        case('ckapppa1')
          description=
     .     'Coating Kappa_perp A1 polynomial factor (W/m/K)'
          nz=1
          call wall_copy(coating_Kappa_perp_A1,1.0_R8,nz,textbl)
        case('ckapppa2')
          description=
     .     'Coating Kappa_perp A2 polynomial factor (W/m/K^{2})'
          nz=1
          call wall_copy(coating_Kappa_perp_A2,1.0_R8,nz,textbl)
        case('ckapppa3')
          description=
     .     'Coating Kappa_perp A3 polynomial factor (W/m/K^{3})'
          nz=1
          call wall_copy(coating_Kappa_perp_A3,1.0_R8,nz,textbl)
        case('ckapppa4')
          description=
     .     'Coating Kappa_perp A4 polynomial factor (W/m/K^{4})'
          nz=1
          call wall_copy(coating_Kappa_perp_A4,1.0_R8,nz,textbl)
        case('ctgdns')
          description='Coating mass density (kg/m^{3})'
          nz=1
          call wall_copy(coating_target_density,1.0_R8,nz,textbl)
        case('c_cp_b1')
          description='Coating heat capacity B1 factor (J/g/K)'
          nz=1
          call wall_copy(coating_Cp_B1,1.0_R8,nz,textbl)
        case('c_cp_b2')
          description='Coating heat capacity B2 factor (J/g/K^{2})'
          nz=1
          call wall_copy(coating_Cp_B2,1.0_R8,nz,textbl)
        case('c_cp_b3')
          description='Coating heat capacity B3 factor (J.K/g)'
          nz=1
          call wall_copy(coating_Cp_B3,1.0_R8,nz,textbl)
        case('erel')
          description='Chemical sputtering factor erel'
          nz=1
          call wall_copy(matsurf_chemsput_erel,1.0_R8,nz,textbl)
        case('bbdyemis')
          description='Blackbody emissivity'
          nz=1
          call wall_copy(matsurf_blackbody_emiss,1.0_R8,nz,textbl)
        case('vapprs_a')
          description='Vapor pressure A factor (K)'
          nz=1
          call wall_copy(matsurf_vapprs_A,1.0_R8,nz,textbl)
        case('vapprs_b')
          description='Vapor pressure B factor'
          nz=1
          call wall_copy(matsurf_vapprs_B,1.0_R8,nz,textbl)
        case('vapprs')
          description='Vapor pressure (Pa)'
          nz=1
          do ibnd = 1, nbnd
            data(ibnd)=10.0e0_R8**
     .       (matsurf_vapprs_A(ibnd)/target_temp(ibnd,1) + 
     .        matsurf_vapprs_B(ibnd))
          enddo
          call wall_copy(data,1.0_R8,nz,textbl)
        case('qevap')
          description='Energy lost by vaporisation (J/m^{2}/s)'
          nz=1
          do ibnd = 1, nbnd
            data(ibnd)=qevap(ibnd,ns,target_temp(ibnd,1))
          enddo
          call wall_copy(data,1.0_R8,nz,textbl)
        case('qrad')
          description='Energy lost by blackbody radiation (J/m^{2}/s)'
          nz=1
          do ibnd = 1, nbnd
            data(ibnd)=qrad(ibnd,target_temp(ibnd,1))
          enddo
          call wall_copy(data,1.0_R8,nz,textbl)
        case('stckcoef')
          description='Sticking coefficient'
          nz=1
          call wall_copy(matsurf_stick_coef,1.0_R8,nz,textbl)
        case('sublheat')
          description='Heat of sublimation (J/kg)'
          nz=1
          call wall_copy(matsurf_heat_sublim,1.0_R8,nz,textbl)
        case('lkappa')
          description='Layer heat conductivity (W/m/K)'
          nz=1
          do ibnd = 1, nbnd
            data(ibnd)=layer_Kappa(ibnd,target_temp(ibnd,1))
          end do
          call wall_copy(data,1.0_R8,nz,textbl)
        case('lrho')
          description='Layer mass density (kg/m^{3})'
          nz=1
          do ibnd = 1, nbnd
            data(ibnd)=layer_mass_density(ibnd)
          end do
          call wall_copy(data,1.0_R8,nz,textbl)
        case('lcp')
          description='Layer heat capacity (J/kg/K)'
          nz=1
          do ibnd = 1, nbnd
            data(ibnd)=layer_Cp(ibnd,target_temp(ibnd,1))
          end do
          call wall_copy(data,1.0_R8,nz,textbl)
        case('lemis')
          description='Layer blackbody emissivity'
          nz=1
          do ibnd = 1, nbnd
            data(ibnd)=layer_emissivity(ibnd)
          end do
          call wall_copy(data,1.0_R8,nz,textbl)
        case('lthick')
          description='Layer thickness (m)'
          nz=1
          do ibnd = 1, nbnd
            data(ibnd)=layer_thickness(ibnd)
          end do
          call wall_copy(data,1.0_R8,nz,textbl)
        case('ldepth')
          description='Layer depth (m)'
          nz=1
          do ibnd = 1, nbnd
            data(ibnd)=layer_depth(ibnd)
          end do
          call wall_copy(data,1.0_R8,nz,textbl)
        case('sputfrac')
          description='Layer fraction'
          nz=ntargsp+1
          do ibnd = 1, nbnd
            if (sput_frac_flag.eq.0) then
              wall_tmp(ibnd,0) = 1.0_R8
              wall_tmp(ibnd,1:ntargsp) = 0.0_R8
            else
c the base material ("0") is assumed to have an unnormalised weight 
c sput_mixed_alpha (default 1)
              wall_tmp(ibnd,0)=sput_mixed_alpha
c the remaining constituents each have an unnormalised weight equal
c to the number of monolayers deposited raised to the power 
c sput_mixed_beta (default 1)
              do itrack = 1, ntargsp
                wall_tmp(ibnd,0)=wall_tmp(ibnd,0)+
     1           monolayer_deposition(ibnd,itrack)**sput_mixed_beta
              end do
c the remaining constituents are then normalised by the sum of weights
              do itrack=1,ntargsp
                wall_tmp(ibnd,itrack)=
     1           monolayer_deposition(ibnd,itrack)**sput_mixed_beta/
     &           wall_tmp(ibnd,0)
              enddo
c the base material weight is calculated to ensure that the sum of 
c normalised weights is 1
              wall_tmp(ibnd,0)=1.0_R8
              do itrack=1,ntargsp
                wall_tmp(ibnd,0)=wall_tmp(ibnd,0)-wall_tmp(ibnd,itrack)
              enddo
            endif
          enddo
          call wall_copy(wall_tmp,1.0_R8,nz,texttrack)
        case default WVALUE
          cflag=.false.
      END SELECT WVALUE

      if(cflag) then
        call customise_label(description,wall_label(iwwrk))
        call set_wallbarlabel(wall_label(iwwrk))
        if (b2plot_debug.ge.1)
     .   write(*,*) 'Copied ',trim(cval),' to wall stack.'
      endif
      return

      end subroutine select_wall_stack

      subroutine select_wall_command

      cflag=.true.
      WFVALUE: select case (cval)
        case('wall')    ! adds wall stack data to contour plot
          plotting_wall=.true.
          call set_wall_on(values(14),values(15),wlinlog)
        case('wall_ank') ! add Kukushkin wall data to contour plot
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          plotting_wall=istack(iistack) .ne. 0
          call wallon(.true.,istack(iistack),values(14),values(15),
     .     wlinlog)
          call pop_stack(stack_req)
        case('logw')    ! log wall elements
          if (mod(wlinlog,2).eq.1) wlinlog=wlinlog+1
          call set_wlinlog(wlinlog)
        case('linw')    ! linear wall elements
          if (wlinlog.eq.2.or.wlinlog.eq.4) wlinlog=wlinlog-1
          call set_wlinlog(wlinlog)
        case('logxw')   ! log x wall elements plot
          if (wlinlog.le.2) wlinlog=wlinlog+2
          call set_wlinlog(wlinlog)
        case('linxw')   ! linear x wall elements plot
          if (wlinlog.gt.2) wlinlog=wlinlog-2
          call set_wlinlog(wlinlog)
        case('plwl')    ! plot in wall domain
          call walloff
          call plot_wall_elements
        case('sumw')    ! sum over walls
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          call sumw()
          graphlabel(iwwrk)=graphlabel(iwwrk)
     1                 (1:lnblnk(graphlabel(iwwrk)))//' summed '
          call set_wallbarlabel(wall_label(iwwrk))
        case('intw')    ! integrate over walls
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          call intw()
          graphlabel(iwwrk)=graphlabel(iwwrk)
     1                 (1:lnblnk(graphlabel(iwwrk)))//' integrated '
          call set_wallbarlabel(wall_label(iwwrk))
        case('intw*')    ! integrate over walls (no account of segments)
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          call intwstar()
          graphlabel(iwwrk)=graphlabel(iwwrk)
     1                 (1:lnblnk(graphlabel(iwwrk)))//' integrated '
          call set_wallbarlabel(wall_label(iwwrk))
        case('pswa')     ! per wall area
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          call pswa
        case('wclvs')    ! Customise Wall Contour Levels
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          stack_req=(/1,istack(iistack),0,0,0,0/)
          if(bad_args(stack_req)) return
          call wclvs(istack(iistack),rstack,irstack)
          do iclv=1,istack(iistack)
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('wclvd')    ! Customise Wall Contour Level Delta
          stack_req=(/0,3,0,0,0,0/)
          if(bad_args(stack_req)) return
          call 
     .     wclvd(rstack(irstack-2),rstack(irstack-1),rstack(irstack))
          do iclv=1,3
            dummy=rpop(rstack,irstack)
          end do
        case('wclvn')    ! Customise Wall Contour Level Range/Num
          stack_req=(/1,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          call wclvn(rstack(irstack-1),rstack(irstack),istack(iistack))
          do iclv=1,2
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('wclvp')    ! Customise Wall Contour Level by Proportion
          stack_req=(/2,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          ix=ipop(istack,iistack)
          iy=ipop(istack,iistack)
          call wclvp(ix,iy)
          write(stmp,'(a,2x,a1,2i4,a2)') 
     -     wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk))),'{',ix,iy,' }'
          wall_label(iwwrk)=stmp
          call set_wallbarlabel(wall_label(iwwrk))
        case('wwid')    ! Wall Line Width Factor
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          call walllinewidth(rstack(irstack))
          wwidth = rstack(irstack)
          dummy=rpop(rstack,irstack)
        case('wabs')     ! abs() on wall stack
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do i=-1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=abs(wall_stack(i,is,iwwrk))
            end do
          end do
          wall_label(iwwrk)=
     .      graphlabel(iwwrk)(1:lnblnk(graphlabel(iwwrk)))//' (ABS)'
          call set_wallbarlabel(wall_label(iwwrk))
        case('ddw')     ! d/dw derivative along wall length
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do iseg=1,nsegments
              do i=wall_begin(iseg),wall_end(iseg)
                if(i.eq.wall_begin(iseg).and.i.lt.wall_end(iseg)) then
                  wall_tmp(i,is)=
     1             (wall_stack(i+1,is,iwwrk)-wall_stack(i,is,iwwrk))/
     2             (0.5_R8*(
     .              sqrt((r2(i+1)-r1(i+1))**2+(z2(i+1)-z1(i+1))**2)+
     .              sqrt((r2(i)-r1(i))**2+(z2(i)-z1(i))**2)))
                elseif(i.eq.wall_end(iseg).and.
     &                 i.gt.wall_begin(iseg)) then
                  wall_tmp(i,is)=
     1             (wall_stack(i,is,iwwrk)-wall_stack(i-1,is,iwwrk))/
     2             (0.5_R8*(
     .              sqrt((r2(i-1)-r1(i-1))**2+(z2(i-1)-z1(i-1))**2)+
     .              sqrt((r2(i)-r1(i))**2+(z2(i)-z1(i))**2)))
                elseif(i.ne.wall_end(iseg).and.
     &                 i.ne.wall_begin(iseg)) then
                  wall_tmp(i,is)=
     1             (wall_stack(i+1,is,iwwrk)-wall_stack(i-1,is,iwwrk))/
     .             (0.5_R8*(
     .              sqrt((r2(i-1)-r1(i-1))**2+(z2(i-1)-z1(i-1))**2)+
     .              sqrt((r2(i+1)-r1(i+1))**2+(z2(i+1)-z1(i+1))**2))+
     .              sqrt((r1(i)-r2(i))**2+(z1(i)-z2(i))**2))
                else
                  wall_tmp(i,is) = 0.0_R8
                endif
              end do
            end do
          end do
          wall_stack(1:nbnd,0:nwdim(2,iwwrk)-1,iwwrk)=
     .     wall_tmp(1:nbnd,0:nwdim(2,iwwrk)-1)
          wall_label(iwwrk)=
     .     wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk)))//'/(m) (d/dw)'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wsumz')    ! sum over ns
          stack_req=(/2,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          izstrt=istack(iistack-1)
          izend=istack(iistack)
          if(izstrt .eq. 0 .and. izend .eq. 0) then
            izstrt=0
            izend=nwdim(2,iwwrk)-1
          end if
          if(izend .eq. -1) izend=nwdim(2,iwwrk)-1
          if(izstrt .lt. 0 .or. izend .gt. nwdim(2,iwwrk)-1 .or.
     1       izstrt .gt. izend) then
            write(*,*) 'Parameters supplied to wsumz were invalid - ',
     1                 'valid parameters are:'
            write(*,*) ' 0 <= ',izstrt,' <= ',izend,' <= ',
     1           nwdim(2,iwwrk)-1
            return
          end if
          do is=izstrt+1,izend
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,izstrt,iwwrk)=
     1          wall_stack(i,izstrt,iwwrk)
     2         +wall_stack(i,is,iwwrk)
            end do
          end do
          do i=1,nwdim(1,iwwrk)
            wall_stack(i,0,iwwrk)=
     1       wall_stack(i,izstrt,iwwrk)
          end do
          nwdim(2,iwwrk)=1
          if(izstrt .ne. izend) then
            wall_label(iwwrk)=
     1       wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk))) //' SUM'//
     2       textwwk(izstrt,iwwrk)(1:lnblnk(textwwk(izstrt,iwwrk)))
     3       //' to'//
     4       textwwk(izend,iwwrk)(1:lnblnk(textwwk(izend,iwwrk)))
          else
            wall_label(iwwrk)=
     1       wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk))) //
     2       textwwk(izstrt,iwwrk)(1:lnblnk(textwwk(izstrt,iwwrk)))
          end if
          textwwk(0,iwwrk)=textbl(0)
          call pop_stack(stack_req)
          call set_wallbarlabel(wall_label(iwwrk))
        case('pwne')     ! per ne
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          call pwne
        case('wm+')      ! wall stack addition
          stack_req=(/0,0,0,0,2,0/)
          if(bad_args(stack_req)) return   ! 5900
          if(nwdim(2,iwwrk-1) .ne. nwdim(2,iwwrk) .and.
     1       nwdim(2,iwwrk-1) .ne. 1 .and. nwdim(2,iwwrk) .ne. 1) then
            write(*,*) 'Mismatch in nwdim(2)''s in operands'
            goto 5900
          end if
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk+1)=
     1         wall_stack(i,min(is,nwdim(2,iwwrk-1)-1),iwwrk-1)+
     2         wall_stack(i,min(is,nwdim(2,iwwrk)-1),iwwrk)
            end do
          end do
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk-1)=
     1         wall_stack(i,is,iwwrk+1)
            end do
            if(nwdim(2,iwwrk-1) .gt. nwdim(2,iwwrk)) then
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk-1)
            else
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk)
            end if
          end do
          wall_label(iwwrk-1)='Special'
          nwdim(2,iwwrk-1)=max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))
          iwwrk=iwwrk-1
          call set_wallbarlabel(wall_label(iwwrk))
        case('wm-')      ! wall stack subtraction
          stack_req=(/0,0,0,0,2,0/)
          if(bad_args(stack_req)) return   ! 5900
          if(nwdim(2,iwwrk-1) .ne. nwdim(2,iwwrk) .and.
     1       nwdim(2,iwwrk-1) .ne. 1 .and. nwdim(2,iwwrk) .ne. 1) then
            write(*,*) 'Mismatch in nwdim(2)''s in operands'
            goto 5900
          end if
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk+1)=
     1         wall_stack(i,min(is,nwdim(2,iwwrk-1)-1),iwwrk-1)-
     2         wall_stack(i,min(is,nwdim(2,iwwrk)-1),iwwrk)
            end do
          end do
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk-1)=wall_stack(i,is,iwwrk+1)
            end do
            if(nwdim(2,iwwrk-1) .gt. nwdim(2,iwwrk)) then
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk-1)
            else
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk)
            end if
          end do
          wall_label(iwwrk-1)='Special'
          nwdim(2,iwwrk-1)=max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))
          iwwrk=iwwrk-1
          call set_wallbarlabel(wall_label(iwwrk))
        case('wm*')      ! wall stack multiplication
          stack_req=(/0,0,0,0,2,0/)
          if(bad_args(stack_req)) return   ! 5900
          if(nwdim(2,iwwrk-1) .ne. nwdim(2,iwwrk) .and.
     1       nwdim(2,iwwrk-1) .ne. 1 .and. nwdim(2,iwwrk) .ne. 1) then
            write(*,*) 'Mismatch in nwdim(2)''s in operands'
            goto 5900
          end if
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk+1)=
     1         wall_stack(i,min(is,nwdim(2,iwwrk-1)-1),iwwrk-1)*
     2         wall_stack(i,min(is,nwdim(2,iwwrk)-1),iwwrk)
            end do
          end do
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i=1,nwdim(1,iwwrk)
               wall_stack(i,is,iwwrk-1)=wall_stack(i,is,iwwrk+1)
            end do
            if(nwdim(2,iwwrk-1) .gt. nwdim(2,iwwrk)) then
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk-1)
            else
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk)
            end if
          end do
          wall_label(iwwrk-1)='Special'
          nwdim(2,iwwrk-1)=max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))
          iwwrk=iwwrk-1
          call set_wallbarlabel(wall_label(iwwrk))
        case('wm/')      ! wall stack division
          stack_req=(/0,0,0,0,2,0/)
          if(bad_args(stack_req)) return   ! 5900
          if(nwdim(2,iwwrk-1) .ne. nwdim(2,iwwrk) .and.
     1       nwdim(2,iwwrk-1) .ne. 1 .and. nwdim(2,iwwrk) .ne. 1) then
            write(*,*) 'Mismatch in nwdim(2)''s in operands'
            goto 5900
          end if
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i=1,nwdim(1,iwwrk)
              if(wall_stack(i,min(is,nwdim(2,iwwrk)-1),iwwrk) 
     .         .ne. 0.0_R8) then
               wall_stack(i,is,iwwrk+1)=
     1          wall_stack(i,min(is,nwdim(2,iwwrk-1)-1),iwwrk-1)/
     2          wall_stack(i,min(is,nwdim(2,iwwrk)-1),iwwrk)
              else
                wall_stack(i,is,iwwrk+1)=0.0_R8
              end if
            end do
          end do
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk-1)=wall_stack(i,is,iwwrk+1)
            end do
            if(nwdim(2,iwwrk-1) .gt. nwdim(2,iwwrk)) then
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk-1)
            else
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk)
            end if
          end do
          wall_label(iwwrk-1)='Special'
          nwdim(2,iwwrk-1)=max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))
          iwwrk=iwwrk-1
          call set_wallbarlabel(wall_label(iwwrk))
        case('rwm+')     ! wall stack + real
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return   ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=
     1         wall_stack(i,is,iwwrk)+rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('rwm-')     ! wall stack - real
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return   ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=
     1         wall_stack(i,is,iwwrk)-rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('rwm*')     ! wall stack * real
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return   ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=
     1         wall_stack(i,is,iwwrk)*rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('rwm/')     ! wall stack / real
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return   ! 5900
          if(rstack(irstack).eq.0.0_R8) then
            write (*,*) 'Divide by zero, You lose'
            goto 5900
          endif
          do is=0,nwdim(2,iwwrk)-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=
     1         wall_stack(i,is,iwwrk)/rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('rwm**')    ! wall stack ** real
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return   ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i=1,nwdim(1,iwwrk)
              if (rstack(irstack).ne.0.0) then
                wall_stack(i,is,iwwrk)=
     1           wall_stack(i,is,iwwrk)**rstack(irstack)
              else
                wall_stack(i,is,iwwrk)=1.0_R8
              end if
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wdrop')    ! drop top of wall stack
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return   ! 5900
          iwwrk=iwwrk-1
        case('wdup')     ! duplicate top of wall stack
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return   ! 5900
          if(iwwrk .ge. wklng) call dropwallstackbottom()
          iwwrk=iwwrk+1
          nwdim(1,iwwrk)=nwdim(1,iwwrk-1)
          nwdim(2,iwwrk)=nwdim(2,iwwrk-1)
          do is=0,nwdim(2,iwwrk-1)-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=wall_stack(i,is,iwwrk-1)
            end do
            textwwk(is,iwwrk)=textwwk(is,iwwrk-1)
          end do
          wall_label(iwwrk)=wall_label(iwwrk-1)
          call set_wallbarlabel(wall_label(iwwrk))
        case('wswap')    ! swap top two elements of wall stack
          stack_req=(/0,0,0,0,2,0/)
          if(bad_args(stack_req)) return   ! 5900
          call iswap(nwdim(1,iwwrk),nwdim(1,iwwrk-1))
          call iswap(nwdim(2,iwwrk),nwdim(2,iwwrk-1))
          do is=0,nwdim(2,iwwrk-1)-1
            do i=1,nwdim(1,iwwrk)
              call rswap(wall_stack(i,is,iwwrk),
     1                   wall_stack(i,is,iwwrk-1))
            end do
            call strswap(textwwk(is,iwwrk),textwwk(is,iwwrk-1))
          end do
          call strswap(wall_label(iwwrk),wall_label(iwwrk-1))
          call set_wallbarlabel(wall_label(iwwrk))
        case('wover')    ! w2 w1 --> w2 w1 w2
          stack_req=(/0,0,0,0,2,0/)
          if(bad_args(stack_req)) return ! 5900
          if(iwwrk .ge. wklng) call dropwallstackbottom()
          iwwrk=iwwrk+1
          nwdim(1,iwwrk)=nwdim(1,iwwrk-2)
          nwdim(2,iwwrk)=nwdim(2,iwwrk-2)
          do is=0,nwdim(2,iwwrk-2)-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=wall_stack(i,is,iwwrk-2)
            end do
            textwwk(is,iwwrk)=textwwk(is,iwwrk-2)
          end do
          wall_label(iwwrk)=wall_label(iwwrk-2)
          call set_wallbarlabel(wall_label(iwwrk))
        case('wlist')    ! list contents of wall stack
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do i = iwwrk, 1, -1
            write(*,*) trim(wall_label(i)), 
     &       (trim(textwwk(is,i)),is=0,nwdim(2,i)-1)
          end do
        case('wsel')    ! select wall element
          stack_req=(/1,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          i1=istack(iistack)
          if(1.gt.i1.or.i1.gt.nwdim(1,iwwrk)) then
            write(*,*) 'Faulty parameter i1 = ',i1,' for wsel!'
            return
          endif
          do is=0,nwdim(2,iwwrk)-1
            do i=1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=wall_stack(i1,is,iwwrk)
            end do
          end do
          call pop_stack(stack_req)
          call set_wallbarlabel(wall_label(iwwrk))
        case('wzsel')    ! select species in wall stack
          stack_req=(/1,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          is1=istack(iistack)
          if(0.gt.is1.or.is1.ge.nwdim(2,iwwrk)) then
            write(*,*) 'Faulty parameter is = ',is1,' for wzsel!'
            return
          endif
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=wall_stack(i,is1,iwwrk)
            end do
            textwwk(is,iwwrk)=textwwk(is1,iwwrk)
          end do
          nwdim(2,iwwrk)=1
          call pop_stack(stack_req)
        case('wmmin')    ! min(w1,w2)
          stack_req=(/0,0,0,0,2,0/) 
          if(bad_args(stack_req)) return ! 5900
          if(nwdim(2,iwwrk-1) .ne. nwdim(2,iwwrk) .and.
     1       nwdim(2,iwwrk-1) .ne. 1 .and. nwdim(2,iwwrk) .ne. 1) then
            write(*,*) 'Mismatch in nwdim(2)''s in operands'
            goto 5900
          end if
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk+1)=
     1         min(wall_stack(i,min(is,nwdim(2,iwwrk-1)-1),iwwrk-1),
     2             wall_stack(i,min(is,nwdim(2,iwwrk)-1),iwwrk))
            end do
          end do
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk-1)=wall_stack(i,is,iwwrk+1)
            end do
            if(nwdim(2,iwwrk-1) .gt. nwdim(2,iwwrk)) then
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk-1)
            else
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk)
            end if
          end do
          wall_label(iwwrk-1)='Special'
          nwdim(2,iwwrk-1)=max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))
          iwwrk=iwwrk-1
          call set_wallbarlabel(wall_label(iwwrk))
        case('wmmax')    ! max(w1,w2)
          stack_req=(/0,0,0,0,2,0/)
          if(bad_args(stack_req)) return ! 5900
          if(nwdim(2,iwwrk-1) .ne. nwdim(2,iwwrk) .and.
     1       nwdim(2,iwwrk-1) .ne. 1 .and. nwdim(2,iwwrk) .ne. 1) then
            write(*,*) 'Mismatch in nwdim(2)''s in operands'
            goto 5900
          end if
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk+1)=
     1         max(wall_stack(i,min(is,nwdim(2,iwwrk-1)-1),iwwrk-1),
     2             wall_stack(i,min(is,nwdim(2,iwwrk)-1),iwwrk))
            end do
          end do
          do is=0,max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk-1)=wall_stack(i,is,iwwrk+1)
            end do
            if(nwdim(2,iwwrk-1) .gt. nwdim(2,iwwrk)) then
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk-1)
            else
              textwwk(is,iwwrk-1)=textwwk(is,iwwrk)
            end if
          end do
          wall_label(iwwrk-1)='Special'
          nwdim(2,iwwrk-1)=max(nwdim(2,iwwrk-1),nwdim(2,iwwrk))
          iwwrk=iwwrk-1
          call set_wallbarlabel(wall_label(iwwrk))
        case('rwmmi')    ! min(r1,w1)
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=
     1         min(wall_stack(i,is,iwwrk),rstack(irstack)) 
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('rwmma')    ! max (r1,w1)
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=
     1         max(wall_stack(i,is,iwwrk),rstack(irstack))
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wmlt')    ! w<r
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              if(wall_stack(i,is,iwwrk).lt.
     1                 rstack(irstack)) then
                wall_stack(i,is,iwwrk)=1.0_R8
              else
                wall_stack(i,is,iwwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wmle')    ! w<=r
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              if(wall_stack(i,is,iwwrk).le.
     1                 rstack(irstack)) then
                wall_stack(i,is,iwwrk)=1.0_R8
              else
                wall_stack(i,is,iwwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wmeq')    ! w=r
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              if(wall_stack(i,is,iwwrk).eq.
     1                 rstack(irstack)) then
                wall_stack(i,is,iwwrk)=1.0_R8
              else
                wall_stack(i,is,iwwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wmne')    ! w!=r
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              if(wall_stack(i,is,iwwrk).ne.
     1                 rstack(irstack)) then
                wall_stack(i,is,iwwrk)=1.0_R8
              else
                wall_stack(i,is,iwwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wmgt')    ! w>r
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              if(wall_stack(i,is,iwwrk).gt.
     1                 rstack(irstack)) then
                wall_stack(i,is,iwwrk)=1.0_R8
              else
                wall_stack(i,is,iwwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wmge')    ! w=>r
          stack_req=(/0,1,0,0,1,0/)
          if(bad_args(stack_req)) return ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              if(wall_stack(i,is,iwwrk).ge.
     1                 rstack(irstack)) then
                wall_stack(i,is,iwwrk)=1.0_R8
              else
                wall_stack(i,is,iwwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wlog')     ! log10(wall)
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              if(wall_stack(i,is,iwwrk).gt.0.0_R8) then
                wall_stack(i,is,iwwrk)=log10(wall_stack(i,is,iwwrk))
              else
                write(*,*) 'Negative argument for ',
     1            i,is,iwwrk,' : Nothing done.'
              end if
            end do
          end do
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wln')      ! ln(wall)
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              if(wall_stack(i,is,iwwrk).gt.0.0_R8) then
                wall_stack(i,is,iwwrk)=log(wall_stack(i,is,iwwrk))
              else
                write(*,*) 'Negative argument for ',
     1           i,is,iwwrk,'. Nothing done.'
              end if
            end do
          end do
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wsqrt')    ! sqrt(wall)
           stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              if(wall_stack(i,is,iwwrk).ge.0.0_R8) then
                wall_stack(i,is,iwwrk)=sqrt(wall_stack(i,is,iwwrk))
              else
                write(*,*) 'Negative argument for ',
     1            i,is,iwwrk,' : Nothing done.'
              end if
            end do
          end do
          call pop_stack(stack_req)
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wexp')     ! exp(wall)
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=exp(wall_stack(i,is,iwwrk))
            end do
          end do
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wzrng')    ! select species range in wall stack
          stack_req=(/2,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          if(0.gt.istack(iistack-1).or.
     &       istack(iistack-1).ge.nwdim(2,iwwrk)) then
            write(*,*) 'Faulty first parameter is = ',istack(iistack-1),
     &       ' for zrng!'
            return
          endif
          if(0.gt.istack(iistack).or.
     &       istack(iistack).ge.nwdim(2,iwwrk)) then
            write(*,*) 'Faulty second parameter is = ',istack(iistack),
     &       ' for zrng!'
            return
          endif
          do is=istack(iistack-1),istack(iistack)
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is-istack(iistack-1),iwwrk)=
     1         wall_stack(i,is,iwwrk)
            end do
            textwwk(is-istack(iistack-1),iwwrk)=textwwk(is,iwwrk)
          end do
          nwdim(2,iwwrk)=istack(iistack)-istack(iistack-1)+1
          call pop_stack(stack_req)
        case('wchs')     ! change sign in wall stack
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=-wall_stack(i,is,iwwrk)
            end do
          end do
          wall_label(iwwrk)=
     1      wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk)))//' (CHS)'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wtos')    ! wall array -> m
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          if(iwrk .eq. wklng) call dropstackbottom()
          iwrk=iwrk+1
          ndim(1,iwrk)=nx
          ndim(2,iwrk)=ny
          ndim(3,iwrk)=nwdim(2,iwwrk)
          work(-1:nx,-1:ny,0:ndim(3,iwrk)-1,iwrk)=0.0_R8
          do i=1,nbnd
            work(imapx(i),imapy(i),0:ndim(3,iwrk)-1,iwrk)=
     .       wall_stack(i,0:ndim(3,iwrk)-1,iwwrk)
          enddo
          values(9)=minval(work(-1:ndim(1,iwrk),-1:ndim(2,iwrk),
     1         0:ndim(3,iwrk)-1,iwrk))
          values(10)=maxval(work(-1:ndim(1,iwrk),-1:ndim(2,iwrk),
     1         0:ndim(3,iwrk)-1,iwrk))          
          textwk(0:ndim(3,iwrk)-1,iwrk)=textwwk(0:ndim(3,iwrk)-1,iwwrk)
          graphlabel(iwrk)=wall_label(iwwrk)
          goto 5059
        case('wtot')    ! wall array --> target
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          if(itwrk .eq. wklng) call droptargetstackbottom()
          itwrk=itwrk+1
          if(nwdim(2,iwwrk).gt.1) 
     >     write(*,*) 'Only the first species of the wall stack ',
     .                'is copied onto the target stack !'
          do i = 1, ndepth_nml
            target_stack(1:nbnd,i,itwrk)=wall_stack(1:nbnd,0,iwwrk)
          enddo
          values(32)=minval(target_stack(1:nbnd,1:ndepth_nml,itwrk))
          values(33)=maxval(target_stack(1:nbnd,1:ndepth_nml,itwrk))          
          target_label(itwrk)=wall_label(iwwrk)
          call set_targetbarlabel(target_label(itwrk))
          goto 5059
        case('wpos')     ! wall stack with 1 if element > 0, else 0
          stack_req=(/0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do i=1,nwdim(1,iwwrk)
              if(wall_stack(i,is,iwwrk) .gt. 0.0_R8) then
                wall_stack(i,is,iwwrk)=1.0_R8
              else
                wall_stack(i,is,iwwrk)=0.0_R8
              endif
            end do
          end do
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wsaw')     ! reads external data onto the wall stack
          stack_req=(/1,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          arg=spop(sstack,isstack)
          ndata=ipop(istack,iistack)
          nz=1
          open(99,file=arg,status='old')
          do i=1,ndata
            read(99,*) r1(i),r2(i),z1(i),z2(i),data(i)
          enddo
          values(5)=min(minval(r1),minval(r2))             	 !rmin
          values(6)=max(maxval(r1),maxval(r2))             	 !rmax
          values(7)=min(minval(z1),minval(z2))             	 !zmin
          values(8)=max(maxval(z1),maxval(z2))             	 !zmax
          write(*,*) 'values',(i,values(i),i=5,8)
          call wall_copy(data,1.0_R8,nz,textbl)
          wall_label(iwwrk)=arg
          call set_wallbarlabel(wall_label(iwwrk))
        case('wallf')    ! add experimental wall data to plot
          stack_req=(/0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          wallfile=spop(sstack,isstack)
          call plot_wallfile
        case('outpw')    ! prints the value of a specific wall element
          stack_req=(/1,0,1,0,1,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          i=ipop(istack,iistack)
          do is=0,nwdim(2,iwwrk)-1
            if(nwdim(2,iwwrk) .eq. 1) then
              write(*,*) stmp(1:lnblnk(stmp)),
     1                 wall_stack(i,is,iwwrk)
              if(writing_to_file) then
                write(nwtf,*) stmp(1:lnblnk(stmp)),                   !srv 21.09.11
     1                   wall_stack(i,is,iwwrk)
              endif
            else
              write(*,*) stmp(1:lnblnk(stmp)),is,
     1                wall_stack(i,is,iwwrk)
              if(writing_to_file) then
                write(nwtf,*) stmp(1:lnblnk(stmp)),is,                   !srv 21.09.11
     1                    wall_stack(i,is,iwwrk)
              endif
            endif
          enddo
        case('wprt')    ! output the top of wall stack to a file
          stack_req=(/0,0,1,0,1,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          open(99,file=stmp)
          do is=0,nwdim(2,iwwrk)-1
            write(99,'(''#'',1x,a,1x,a)')
     1       wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk))),
     2       textwwk(is,iwwrk)(1:lnblnk(textwwk(is,iwwrk)))
            do i=1,nbnd
              write(99,*) i,is,wall_stack(i,is,iwwrk)
            enddo
          enddo
          close(99)
        case('wtab')    ! output a table from the top of wall stack to a file
          stack_req=(/2,0,1,0,1,0/)
          if(bad_args(stack_req)) return
          ix2=ipop(istack,iistack)
          ix1=ipop(istack,iistack)
          if (ix1.eq.0 .and. ix2.eq.0) then
            ix1 = 1
            ix2 = nwdim(1,iwwrk)
          end if
          if (ix1 .lt. 1 .or. ix2 .gt. nwdim(1,iwwrk) .or. 
     1        ix1 .gt. ix2 ) then
            write(*,*) 'Parameters supplied to wtab were invalid - ',
     1                 'valid parameters are:'
            write(*,*) ' 1 <= ',ix1,' <= ',ix2,' <= ',nwdim(1,iwwrk)
            return
          end if
          stmp=spop(sstack,isstack)
          nz=nwdim(2,iwwrk)
          open(99,file=stmp)
          write(hlp_frm,'(a,i2,a)') '(a1,1x,a,1x,',nz,'a13)'
          write(*,*) hlp_frm
          write(99,hlp_frm) '#',
     1       wall_label(iwwrk)(1:lnblnk(wall_label(iwwrk))),
     2       (textwwk(is,iwwrk)(1:lnblnk(textwwk(is,iwwrk))),is=0,nz-1)
          write(hlp_frm,'(a,i2,a)') '(i4,1p,',nz,'(1e12.5,1x))'
          write(*,*) hlp_frm
          do ix=ix1,ix2
            write(99,hlp_frm) ix,wall_stack(ix,0:nz-1,iwwrk)
          enddo
          close(99)
        case('wlld')    ! ANK wall loading
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call b2pwld(istack(iistack),BoRiS)
          call pop_stack(stack_req)
        case('woff')    ! turn off wall plots
          call walloff
        case default WFVALUE
          cflag=.false.
      END SELECT WFVALUE
c
      goto 5900

 5059 continue
      stack_req=(/0,0,0,0,1,0/)
      if(bad_args(stack_req)) return
      iwwrk=iwwrk-1

 5900 continue
      if(cflag) then
        if(b2plot_debug.ge.4)
     1   write(*,*) 'Applied operator ',trim(cval),
     2   ' to wall stack array'
        if(iwwrk.gt.0) then
          values(14)=minval(wall_stack(1:nwdim(1,iwwrk),
     1         0:nwdim(2,iwwrk)-1,iwwrk))
          values(15)=maxval(wall_stack(1:nwdim(1,iwwrk),
     1         0:nwdim(2,iwwrk)-1,iwwrk))
        endif
      endif
      return

      end subroutine select_wall_command

      subroutine select_target_stack

      cflag=.true.
      TFVALUE: select case (cval)
        case('logtrg')    ! log target elements
          if (mod(trglinlog,2).eq.1) trglinlog=trglinlog+1
          call set_trglinlog(trglinlog)
        case('lintrg')    ! linear target elements
          if (trglinlog.eq.2.or.trglinlog.eq.4) trglinlog=trglinlog-1
          call set_trglinlog(trglinlog)
        case('logxtrg')   ! log x target elements plot
          if (trglinlog.le.2) trglinlog=trglinlog+2
          call set_trglinlog(trglinlog)
        case('linxtrg')   ! linear x target elements plot
          if (trglinlog.gt.2) trglinlog=trglinlog-2
          call set_trglinlog(trglinlog)
        case('tclvs')    ! Customize Target Contour Levels
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          stack_req=(/1,istack(iistack),0,0,0,0/)
          if(bad_args(stack_req)) return
          call tclvs(istack(iistack),rstack,irstack)
          do iclv=1,istack(iistack)
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('tclvd')    ! Customise Target Contour Level Delta
          stack_req=(/0,3,0,0,0,0/)
          if(bad_args(stack_req)) return
          call tclvd
     .     (rstack(irstack-2),rstack(irstack-1),rstack(irstack))
          do iclv=1,3
            dummy=rpop(rstack,irstack)
          end do
        case('tclvn')    ! Customise Target Contour Level Range/Num
          stack_req=(/1,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          call tclvn(rstack(irstack-1),rstack(irstack),istack(iistack))
          do iclv=1,2
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('tclvp')    ! Customise Target Contour Level by Proportion
          stack_req=(/2,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          ix=ipop(istack,iistack)
          iy=ipop(istack,iistack)
          call tclvp(ix,iy)
          write(stmp,'(a,2x,a1,2i4,a2)') 
     -     target_label(itwrk)(1:lnblnk(target_label(itwrk))),
     +     '{',ix,iy,' }'
          target_label(itwrk)=stmp
          call set_targetbarlabel(target_label(itwrk))
        case('trg.w')     ! line plot versus target wall
          stack_req=(/1,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          if(1.le.istack(iistack).and.istack(iistack).le.ndepth_nml) 
     >     then
            call trgw(istack(iistack))
          else
            write(*,*) 
     .       'Faulty i = ',istack(iistack),' chosen for trg.w!'
          endif
          call pop_stack(stack_req)
        case('trg.z')     ! line plot versus target depth
          stack_req=(/1,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          if(1.le.istack(iistack).and.istack(iistack).le.nbnd) then
            call trgz(istack(iistack))
          else
            write(*,*) 
     .       'Faulty ibnd = ',istack(iistack),' chosen for trg.z!'
          endif
          call pop_stack(stack_req)
        case('sumtrgw')    ! sum over targets along walls
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call sumtrgw()
          target_label(itwrk)=target_label(itwrk)
     1     (1:lnblnk(target_label(itwrk)))//' summed over w'
          call set_targetbarlabel(target_label(itwrk))
        case('sumtrgz')    ! sum over targets along z
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call sumtrgz()
          target_label(itwrk)=target_label(itwrk)
     1     (1:lnblnk(target_label(itwrk)))//' summed over z'
          call set_targetbarlabel(target_label(itwrk))
        case('inttrgw')    ! integrate over targets along walls
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call inttrgw()
          target_label(itwrk)=target_label(itwrk)
     1     (1:lnblnk(target_label(itwrk)))//' integrated along w'
          call set_targetbarlabel(target_label(itwrk))
        case('inttrgw*')    ! integrate over targets along walls (no account of segments)
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call inttrgwstar()
          target_label(itwrk)=target_label(itwrk)
     1     (1:lnblnk(target_label(itwrk)))//' integrated along w'
          call set_targetbarlabel(target_label(itwrk))
        case('inttrgz')    ! integrate over targets along z
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call inttrgz()
          target_label(itwrk)=target_label(itwrk)
     1     (1:lnblnk(target_label(itwrk)))//' integrated along z'
          call set_targetbarlabel(target_label(itwrk))
        case('trgabs')     ! absolute value of target data
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          target_stack(1:nbnd,:,itwrk)=abs(target_stack(1:nbnd,:,itwrk))
          target_label(itwrk)=target_label(itwrk)
     1     (1:lnblnk(target_label(itwrk)))//' (ABS) '
          call set_targetbarlabel(target_label(itwrk))
        case('trgddw')     ! target d/dw derivative along wall length
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do is=1,ndepth_nml
            do iseg=1,nsegments
              do i=wall_begin(iseg),wall_end(iseg)
                if(i.eq.wall_begin(iseg).and.i.lt.wall_end(iseg)) then
                  trg_tmp(i,is)=
     1        (target_stack(i+1,is,itwrk)-target_stack(i,is,itwrk))/
     2             (0.5_R8*(
     .              sqrt((r2(i+1)-r1(i+1))**2+(z2(i+1)-z1(i+1))**2)+
     .              sqrt((r2(i)-r1(i))**2+(z2(i)-z1(i))**2)))
                elseif(i.eq.wall_end(iseg).and.
     &                 i.gt.wall_begin(iseg)) then
                  trg_tmp(i,is)=
     1        (target_stack(i,is,itwrk)-target_stack(i-1,is,itwrk))/
     2             (0.5_R8*(
     .              sqrt((r2(i-1)-r1(i-1))**2+(z2(i-1)-z1(i-1))**2)+
     .              sqrt((r2(i)-r1(i))**2+(z2(i)-z1(i))**2)))
                elseif(i.ne.wall_end(iseg).and.
     &                 i.ne.wall_begin(iseg)) then
                  trg_tmp(i,is)=
     1        (target_stack(i+1,is,itwrk)-target_stack(i-1,is,itwrk))/
     .             (0.5_R8*(
     .              sqrt((r2(i-1)-r1(i-1))**2+(z2(i-1)-z1(i-1))**2)+
     .              sqrt((r2(i+1)-r1(i+1))**2+(z2(i+1)-z1(i+1))**2))+
     .              sqrt((r1(i)-r2(i))**2+(z1(i)-z2(i))**2))
                else
                  trg_tmp(i,is) = 0.0_R8
                endif
              end do
            end do
          end do
          target_stack(1:nbnd,:,itwrk)=trg_tmp(1:nbnd,:)
          target_label(itwrk)=
     .     target_label(itwrk)(1:lnblnk(target_label(itwrk)))//
     +     '(/m) (d/dw)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgddp')     ! target d/dp derivative along perpendicular direction
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do is=1,ndepth_nml
            do iseg=1,nsegments
              do i=wall_begin(iseg),wall_end(iseg)
                if(i.eq.wall_begin(iseg).and.i.lt.wall_end(iseg)) then
                  trg_tmp(i,is)=
     1        (target_stack(i+1,is,itwrk)-target_stack(i,is,itwrk))/
     2             (0.5_R8*(target_width(i,is)+target_width(i+1,is)))
                elseif(i.eq.wall_end(iseg).and.
     &                 i.gt.wall_begin(iseg)) then
                  trg_tmp(i,is)=
     1        (target_stack(i,is,itwrk)-target_stack(i-1,is,itwrk))/
     2             (0.5_R8*(target_width(i-1,is)+target_width(i,is)))
                elseif(i.ne.wall_end(iseg).and.
     &                 i.ne.wall_begin(iseg)) then
                  trg_tmp(i,is)=
     1        (target_stack(i+1,is,itwrk)-target_stack(i-1,is,itwrk))/
     .             (0.5_R8*(target_width(i-1,is)+target_width(i+1,is))+
     2                      target_width(i,is))
                else
                  trg_tmp(i,is) = 0.0_R8
                endif
              end do
            end do
          end do
          target_stack(1:nbnd,:,itwrk)=trg_tmp(1:nbnd,:)
          target_label(itwrk)=
     .     target_label(itwrk)(1:lnblnk(target_label(itwrk))) //
     +     '(/m) (d/dp)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgddz')     ! target d/dz derivative along target depth
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do i=1,nbnd
            do is=1,ndepth_nml
              if(is.eq.1.and.is.lt.ndepth_nml) then
                trg_tmp(i,is)=
     1           (target_stack(i,is+1,itwrk)-target_stack(i,is,itwrk))/
     2           (target_depth(i,is+1)-target_depth(i,is))
              elseif(is.eq.ndepth_nml.and.is.gt.1) then
                trg_tmp(i,is)=
     1           (target_stack(i,is,itwrk)-target_stack(i,is-1,itwrk))/
     2           (target_depth(i,is)-target_depth(i,is-1))
              elseif(is.gt.1.and.is.lt.ndepth_nml) then
                trg_tmp(i,is)=
     1           (target_stack(i,is+1,itwrk)-
     2            target_stack(i,is-1,itwrk))/
     .           (target_depth(i,is+1)-target_depth(i,is-1))
              else
                trg_tmp(i,is) = 0.0_R8
              endif
            end do
          end do
          target_stack(1:nbnd,:,itwrk)=trg_tmp(1:nbnd,:)
          target_label(itwrk)=
     .     target_label(itwrk)(1:lnblnk(target_label(itwrk)))//
     +     '(/m) (d/dz)'
          call set_targetbarlabel(target_label(itwrk))
        case('rtrgm+')     ! target stack + real
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return   ! 6900
          do i=1,ndepth_nml
            do ibnd=1,nbnd
              target_stack(ibnd,i,itwrk)=
     1         target_stack(ibnd,i,itwrk)+rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('rtrgm-')     ! target stack - real
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return   ! 6900
          do i=1,ndepth_nml
            do ibnd=1,nbnd
              target_stack(ibnd,i,itwrk)=
     1         target_stack(ibnd,i,itwrk)-rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('rtrgm*')     ! target stack * real
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return   ! 6900
          do i=1,ndepth_nml
            do ibnd=1,nbnd
              target_stack(ibnd,i,itwrk)=
     1         target_stack(ibnd,i,itwrk)*rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('rtrgm/')     ! target stack / real
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return   ! 6900
          if(rstack(irstack).eq.0.0_R8) then
            write (*,*) 'Divide by zero, You lose'
            goto 6900
          endif
          do i=1,ndepth_nml
            do ibnd=1,nbnd
              target_stack(ibnd,i,itwrk)=
     1         target_stack(ibnd,i,itwrk)/rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('rtrgm**')    ! target stack ** real
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return   ! 6900
          do i=1,ndepth_nml
            do ibnd=1,nbnd
              if (rstack(irstack).ne.0.0_R8) then
                target_stack(ibnd,i,itwrk)=
     1           target_stack(ibnd,i,itwrk)**rstack(irstack)
              else
                target_stack(ibnd,i,itwrk)=1.0_R8
              end if
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trgtemp')    ! target temperature (K)
          call target_copy(target_temp,1.0_R8)
          target_label(itwrk)='Target temperature (K)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgdrop')    ! drop top of target stack
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          itwrk=itwrk-1
        case('trgwsel')    ! select wall element in target stack
          stack_req=(/1,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          i1=istack(iistack)
          if(1.gt.i1.or.i1.gt.nbnd) then
            write(*,*) 'Faulty parameter i1 = ',i1,' for trgwsel!'
            return
          endif
          do i=1,ndepth_nml
            do ibnd=1,nbnd
              target_stack(ibnd,i,itwrk)=target_stack(i1,i,itwrk)
            end do
          end do
          call pop_stack(stack_req)
        case('trgzsel')    ! select depth element in target stack
          stack_req=(/1,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          i1=istack(iistack)
          if(1.gt.i1.or.i1.gt.ndepth_nml) then
            write(*,*) 'Faulty parameter i1 = ',i1,' for trgzsel!'
            return
          endif
          do i=1,ndepth_nml
            do ibnd=1,nbnd
              target_stack(ibnd,i,itwrk)=target_stack(ibnd,i1,itwrk)
            end do
          end do
          call pop_stack(stack_req)
        case('rtrgmmi')    ! min(r1,trg)
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return ! 6900
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              target_stack(ibnd,i,itwrk)=
     1         min(target_stack(ibnd,i,itwrk),rstack(irstack)) 
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('rtrgmma')    ! max (r1,trg)
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return ! 6900
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              target_stack(ibnd,i,itwrk)=
     1         max(target_stack(ibnd,i,itwrk),rstack(irstack))
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trglt')    ! trg<r
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return ! 6900
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              if(target_stack(ibnd,i,itwrk).lt.
     1                 rstack(irstack)) then
                target_stack(ibnd,i,itwrk)=1.0_R8
              else
                target_stack(ibnd,i,itwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trgle')    ! trg<=r
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return ! 6900
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              if(target_stack(ibnd,i,itwrk).le.
     1                 rstack(irstack)) then
                target_stack(ibnd,i,itwrk)=1.0_R8
              else
                target_stack(ibnd,i,itwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trgeq')    ! trg=r
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return ! 6900
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              if(target_stack(ibnd,i,itwrk).eq.
     1                 rstack(irstack)) then
                target_stack(ibnd,i,itwrk)=1.0_R8
              else
                target_stack(ibnd,i,itwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trgne')    ! trg!=r
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return ! 6900
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              if(target_stack(ibnd,i,itwrk).ne.
     1                 rstack(irstack)) then
                target_stack(ibnd,i,itwrk)=1.0_R8
              else
                target_stack(ibnd,i,itwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trggt')    ! trg>r
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return ! 6900
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              if(target_stack(ibnd,i,itwrk).gt.
     1                 rstack(irstack)) then
                target_stack(ibnd,i,itwrk)=1.0_R8
              else
                target_stack(ibnd,i,itwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trgge')    ! trg>=r
          stack_req=(/0,1,0,0,0,1/)
          if(bad_args(stack_req)) return ! 6900
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              if(target_stack(ibnd,i,itwrk).ge.
     1                 rstack(irstack)) then
                target_stack(ibnd,i,itwrk)=1.0_R8
              else
                target_stack(ibnd,i,itwrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trglog')     ! log10(target)
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              if(target_stack(ibnd,i,itwrk).gt.0) then
                target_stack(ibnd,i,itwrk)=
     1           log10(target_stack(ibnd,i,itwrk))
              else
                write(*,*) 'Negative argument for ',
     1            ibnd,i,' : Nothing done.'
              end if
            end do
          end do
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trgln')      ! ln(target)
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              if(target_stack(ibnd,i,itwrk).gt.0.0_R8) then
                target_stack(ibnd,i,itwrk)=
     1           log(target_stack(ibnd,i,itwrk))
              else
                write(*,*) 'Negative argument for ',
     1           ibnd,i,' : Nothing done.'
              end if
            end do
          end do
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trgsqrt')    ! sqrt(target)
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              if(target_stack(ibnd,i,itwrk).ge.0.0_R8) then
                target_stack(ibnd,i,itwrk)=
     1           sqrt(target_stack(ibnd,i,itwrk))
              else
                write(*,*) 'Negative argument for ',
     1            ibnd,i,' : Nothing done.'
              end if
            end do
          end do
          call pop_stack(stack_req)
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trgexp')     ! exp(target)
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              target_stack(ibnd,i,itwrk)=
     1         exp(target_stack(ibnd,i,itwrk))
            end do
          end do
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trgchs')     ! change sign in target stack
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              target_stack(ibnd,i,itwrk)=-target_stack(ibnd,i,itwrk)
            end do
          end do
          target_label(itwrk)=
     1      target_label(itwrk)(1:lnblnk(target_label(itwrk)))//' (CHS)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgtow')    ! target array --> wall
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          if(iwwrk .eq. wklng) call dropwallstackbottom()
          iwwrk=iwwrk+1
          nwdim(1,iwwrk)=nbnd
          nwdim(2,iwwrk)=1
          do ibnd=1,nbnd
            wall_stack(ibnd,0,iwwrk)=target_stack(ibnd,1,itwrk)
          enddo
          values(14)=minval(wall_stack(1:nbnd,0,iwwrk))
          values(15)=maxval(wall_stack(1:nbnd,0,iwwrk))          
          textwwk(0,iwwrk)=textbl(0)
          wall_label(iwwrk)=target_label(itwrk)
          call set_wallbarlabel(target_label(itwrk))
        case('trgpos')     ! target stack with 1 if element > 0, else 0
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do i=1,ndepth_nml
            do ibnd=1,nbnd
              if(target_stack(ibnd,i,itwrk) .gt. 0.0_R8) then
                target_stack(ibnd,i,itwrk)=1.0_R8
              else
                target_stack(ibnd,i,itwrk)=0.0_R8
              endif
            end do
          end do
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('outptrg')    ! prints the value of a specific target element
          stack_req=(/2,0,1,0,0,1/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          i=ipop(istack,iistack)
          ibnd=ipop(istack,iistack)
          write(*,*) stmp(1:lnblnk(stmp)),target_stack(ibnd,i,itwrk)
          if(writing_to_file) then
           write(nwtf,*) stmp(1:lnblnk(stmp)),target_stack(ibnd,i,itwrk) !srv 21.09.11
          endif
        case('trgprt')    ! output the target stack to a file
          stack_req=(/0,0,1,0,0,1/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          open(99,file=stmp)
          write(99,'(''#'',1x,a)') 
     .     target_label(itwrk)(1:lnblnk(target_label(itwrk)))
          do i=1,ndepth_nml
            do ibnd=1,nbnd
              write(99,*) ibnd,i,target_stack(ibnd,i,itwrk)
            enddo
          enddo
          close(99)
        case('trgtab')    ! output a table from the top of target stack to a file
          stack_req=(/4,0,1,0,0,1/)
          if(bad_args(stack_req)) return
          iy2=ipop(istack,iistack)
          iy1=ipop(istack,iistack)
          if (iy1.eq.0 .and. iy2.eq.0) then
            iy1 = 1
            iy2 = ndepth_nml
          end if
          ix2=ipop(istack,iistack)
          ix1=ipop(istack,iistack)
          if (ix1.eq.0 .and. ix2.eq.0) then
            ix1 = 1
            ix2 = nbnd
          end if
          if (ix1 .lt. 1 .or. ix2 .gt. nbnd .or.
     &        iy1 .lt. 1 .or. iy2 .gt. ndepth_nml .or.
     1        ix1 .gt. ix2 .or. iy1 .gt. iy2) then
            write(*,*) 'Parameters supplied to trgtab were invalid - ',
     1                 'valid parameters are:'
            write(*,*) ' 1 <= ',ix1,' <= ',ix2,' <= ',nbnd
            write(*,*) ' 1 <= ',iy1,' <= ',iy2,' <= ',ndepth_nml
            return
          end if
          stmp=spop(sstack,isstack)
          open(99,file=stmp)
          write(99,*) '#',
     1       target_label(itwrk)(1:lnblnk(target_label(itwrk)))
          do iy=iy1,iy2
            do ix=ix1,ix2
              write(99,'(2i4,1pe12.5)') ix,iy,target_stack(ix,iy,itwrk)
            enddo
          enddo
          close(99)
        case('trgl')    ! ANK target loading
          stack_req=(/2,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call b2ptrg(istack,iistack)
          call pop_stack(stack_req)
        case('trgsurf')    ! plot target plot
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call trgsurf(1)
        case('kappa_ab')   ! target heat conductivity computed as 1/(AT+B)
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              if (target_depth(ibnd,i).ge.coating_thickness(i)) then
                trg_tmp(ibnd,i) = 1.0_R8 /
     .           (matbulk_Kappa_KA(ibnd)*target_temp(ibnd,i) +
     .            matbulk_Kappa_KB(ibnd))
              else
                trg_tmp(ibnd,i) = 1.0_R8 /
     .           (coating_Kappa_KA(ibnd)*target_temp(ibnd,i) +
     .            coating_Kappa_KB(ibnd))
              end if
            end do
          end do
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='Kappa from 1/(AT+B) (W/m/K)'
          call set_targetbarlabel(target_label(itwrk))
        case('kappa')   ! target heat conductivity computed from polynomial fit
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              trg_tmp(ibnd,i) = Kappa(ibnd,i,target_temp(ibnd,i))
            end do
          end do
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='Target heat conductivity (W/m/K)'
          call set_targetbarlabel(target_label(itwrk))
        case('dkappadt')   ! target heat conductivity derivative
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              trg_tmp(ibnd,i)=dKappabydT(ibnd,i,target_temp(ibnd,i))
            end do
          end do
          call target_copy(trg_tmp,1.0_R8)
          call customise_label
     .     ('(d Kappa/dT) (W/m/K^{2})',target_label(itwrk))
          call set_targetbarlabel(target_label(itwrk))
        case('kappap')   ! perp. target heat conductivity computed from polynomial fit
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              trg_tmp(ibnd,i) =
     .         Kappa_perp(ibnd,i,target_temp(ibnd,i))
            end do
          end do
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)=
     .     'Perpendicular target heat conductivity (W/m/K)'
          call set_targetbarlabel(target_label(itwrk))
        case('kappapdt')   ! perp. target heat conductivity derivative
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              trg_tmp(ibnd,i)=
     .         dKappa_perpbydT(ibnd,i,target_temp(ibnd,i))
            end do
          end do
          call target_copy(trg_tmp,1.0_R8)
          call customise_label(
     .     '(d Kappa_perp/dT) (W/m/K^{2})',target_label(itwrk))
          call set_targetbarlabel(target_label(itwrk))
        case('cp')   ! target specific heat capacity
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              trg_tmp(ibnd,i)=Cp(ibnd,i,target_temp(ibnd,i))
            end do
          end do
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='Target heat capacity (J/kg/K)'
          call set_targetbarlabel(target_label(itwrk))
        case('rho')   ! target mass density          
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              trg_tmp(ibnd,i)=mass_density(ibnd,i)
            end do
          end do
          call target_copy(trg_tmp,1.0_R8)
          call customise_label
     .     ('Target mass density (kg/m^{3})',target_label(itwrk))
          call set_targetbarlabel(target_label(itwrk))
        case('trgwidth')   ! target element width
          call target_copy(target_width,1.0_R8)
          target_label(itwrk)='Target element width (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgdepth')   ! target element depth position
          call target_copy(target_depth,1.0_R8)
          target_label(itwrk)='Target element depth (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgvol')     ! target element volume
          call target_copy(target_vol,1.0_R8)
          call customise_label
     .     ('Target element volume (m^{3})',target_label(itwrk))
          call set_targetbarlabel(target_label(itwrk))
        case('trgdup')     ! duplicate top of target stack
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          if(itwrk .ge. wklng) call droptargetstackbottom()
          itwrk=itwrk+1
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              target_stack(ibnd,i,itwrk)=target_stack(ibnd,i,itwrk-1)
            end do
          end do
          target_label(itwrk)=target_label(itwrk-1)
          call set_targetbarlabel(target_label(itwrk))
        case('trgswap')    ! swap top two elements of target stack
          stack_req=(/0,0,0,0,0,2/)
          if(bad_args(stack_req)) return
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              call rswap(target_stack(ibnd,i,itwrk),
     1                   target_stack(ibnd,i,itwrk-1))
            end do
          end do
          call strswap(target_label(itwrk),target_label(itwrk-1))
          call set_targetbarlabel(target_label(itwrk))
        case('trgover')    ! t2 t1 --> t2 t1 t2
          stack_req=(/0,0,0,0,0,2/)
          if(bad_args(stack_req)) return
          if(iwwrk .ge. wklng) call droptargetstackbottom()
          itwrk=itwrk+1
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              target_stack(ibnd,i,itwrk)=target_stack(ibnd,i,itwrk-2)
            end do
          end do
          target_label(itwrk)=target_label(itwrk-2)
          call set_targetbarlabel(target_label(itwrk))
        case('trglist')    ! list contents of target stack
          stack_req=(/0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do i = itwrk, 1, -1
            write(*,'(a)') trim(target_label(i))
          end do
        case('trgmmin')    ! min(t1,t2)
          stack_req=(/0,0,0,0,0,2/) 
          if(bad_args(stack_req)) return ! 5900
          do ibnd=1,nbnd
            do i = 1, ndepth_nml
              target_stack(ibnd,i,itwrk+1)=
     1         min(target_stack(ibnd,i,itwrk-1),
     2             target_stack(ibnd,i,itwrk))
            end do
          end do
          do ibnd=1,nbnd
            do i = 1, ndepth_nml
              target_stack(ibnd,i,itwrk-1)=target_stack(ibnd,i,itwrk+1)
            end do
          end do
          target_label(itwrk-1)='Special'
          itwrk=itwrk-1
          call set_targetbarlabel(target_label(itwrk))
        case('trgmmax')    ! max(t1,t2)
          stack_req=(/0,0,0,0,0,2/)
          if(bad_args(stack_req)) return ! 5900
          do ibnd=1,nbnd
            do i = 1, ndepth_nml
              target_stack(ibnd,i,itwrk+1)=
     1         max(target_stack(ibnd,i,itwrk-1),
     2             target_stack(ibnd,i,itwrk))
            end do
          end do
          do ibnd=1,nbnd
            do i = 1, ndepth_nml
              target_stack(ibnd,i,itwrk-1)=target_stack(ibnd,i,itwrk+1)
            end do
          end do
          target_label(itwrk-1)='Special'
          itwrk=itwrk-1
          call set_targetbarlabel(target_label(itwrk))
        case('trgr1')   ! Target element R1 position (m)
          call target_copy(targetcr(:,:,1),1.0_R8)
          target_label(itwrk)='Target element R1 position (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgr2')   ! Target element R2 position (m)
          call target_copy(targetcr(:,:,2),1.0_R8)
          target_label(itwrk)='Target element R2 position (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgz1')   ! Target element Z1 position (m)
          call target_copy(targetcz(:,:,1),1.0_R8)
          target_label(itwrk)='Target element Z1 position (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgz2')   ! Target element Z2 position (m)
          call target_copy(targetcz(:,:,2),1.0_R8)
          target_label(itwrk)='Target element Z2 position (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgnorth')  ! North boundary
          do ibnd = 1, nbnd
            if(topiy(imapx(ibnd),imapy(ibnd)).eq.ny+1) then
              trg_tmp(ibnd,1:ndepth_nml) = 1.0_R8
            else
              trg_tmp(ibnd,1:ndepth_nml) = 0.0_R8
            endif
          enddo
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='North boundary'
          call set_targetbarlabel(target_label(itwrk))
        case('trgsouth')  ! South boundary
          do ibnd = 1, nbnd
            if(bottomiy(imapx(ibnd),imapy(ibnd)).eq.-2) then
              trg_tmp(ibnd,1:ndepth_nml) = 1.0_R8
            else
              trg_tmp(ibnd,1:ndepth_nml) = 0.0_R8
            endif
          enddo
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='South boundary'
          call set_targetbarlabel(target_label(itwrk))
        case('trgwest')  ! West boundary
          do ibnd = 1, nbnd
            if(leftix(imapx(ibnd),imapy(ibnd)).eq.-2) then
              trg_tmp(ibnd,1:ndepth_nml) = 1.0_R8
            else
              trg_tmp(ibnd,1:ndepth_nml) = 0.0_R8
            endif
          enddo
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='West boundary'
          call set_targetbarlabel(target_label(itwrk))
        case('trgeast')  ! East boundary
          do ibnd = 1, nbnd
            if(rightix(imapx(ibnd),imapy(ibnd)).eq.nx+1) then
              trg_tmp(ibnd,1:ndepth_nml) = 1.0_R8
            else
              trg_tmp(ibnd,1:ndepth_nml) = 0.0_R8
            endif
          enddo
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='East boundary'
          call set_targetbarlabel(target_label(itwrk))
        case('trgseg')  ! Target segment
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(abs(istack(iistack)).gt.nsegments) then
            write(*,*) 'No such wall segment'
            return
          endif
          if(istack(iistack) .eq. 0) then
            do i = 1, nsegments
              trg_tmp(wall_begin(i):wall_end(i),1:ndepth_nml) = real(i)
            enddo
          elseif(istack(iistack) .gt. 0) then
            trg_tmp(wall_begin(1):
     .              wall_end(istack(iistack)-1),1:ndepth_nml) = 0.0_R8
            trg_tmp(wall_begin(istack(iistack)):
     .              wall_end(istack(iistack)),1:ndepth_nml) = 1.0_R8
            trg_tmp(wall_begin(istack(iistack)+1):
     .              wall_end(nsegments),1:ndepth_nml) = 0.0_R8
          else
            trg_tmp(wall_begin(1):
     .              wall_end(abs(istack(iistack))-1),1:ndepth_nml) = 
     .                   1.0_R8
            trg_tmp(wall_begin(abs(istack(iistack))):
     .              wall_end(abs(istack(iistack))),1:ndepth_nml) = 
     .                   0.0_R8
            trg_tmp(wall_begin(abs(istack(iistack))+1):
     .              wall_end(nsegments),1:ndepth_nml) = 1.0_R8
          endif
          write(*,*) '"trgseg" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='Target segment'
          call set_targetbarlabel(target_label(itwrk))
        case ('trgiz')   ! trgiz
          do is = 1, ndepth_nml
            trg_tmp(1:nbnd,is) = real(is)
          end do
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='Target depth index'
          call set_targetbarlabel(target_label(itwrk))
        case ('itrglt')  ! itrglt
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          is = istack(iistack)
          if(is.gt.nbnd .or. is.le.0) then
            write(*,*) 'No such wall element'
            return
          endif
          trg_tmp(1:is-1,1:ndepth_nml) = 1.0_R8
          trg_tmp(is:nbnd,1:ndepth_nml) = 0.0_R8
          write(*,*) '"itrglt" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='itrglt'
          call set_targetbarlabel(target_label(itwrk))
        case ('itrgeq')  ! itrgeq
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          is = istack(iistack)
          if(is.gt.nbnd .or. is.le.0) then
            write(*,*) 'No such wall element'
            return
          endif
          trg_tmp(1:nbnd,1:ndepth_nml) = 0.0_R8
          trg_tmp(is,1:ndepth_nml) = 1.0_R8
          write(*,*) '"itrgeq" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='itrgeq'
          call set_targetbarlabel(target_label(itwrk))
        case ('itrggt')  ! itrggt
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          is = istack(iistack)
          if(is.gt.nbnd .or. is.le.0) then
            write(*,*) 'No such wall element'
            return
          endif
          trg_tmp(1:is,1:ndepth_nml) = 0.0_R8
          trg_tmp(is+1:nbnd,1:ndepth_nml) = 1.0_R8
          write(*,*) '"itrggt" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='itrggt'
          call set_targetbarlabel(target_label(itwrk))
        case ('ztrglt')  ! ztrglt
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp = istack(iistack)
          if(itmp.gt.ndepth_nml .or. itmp.le.0) then
            write(*,*) 'No such wall depth index'
            return
          endif
          trg_tmp(1:nbnd,1:itmp-1) = 1.0_R8
          trg_tmp(1:nbnd,itmp:ndepth_nml) = 0.0_R8
          write(*,*) '"ztrglt" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='ztrglt'
          call set_targetbarlabel(target_label(itwrk))
        case ('ztrgeq')  ! ztrgeq
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp = istack(iistack)
          if(itmp.gt.ndepth_nml .or. itmp.le.0) then
            write(*,*) 'No such wall depth index'
            return
          endif
          trg_tmp(1:nbnd,1:ndepth_nml) = 0.0_R8
          trg_tmp(1:nbnd,itmp) = 1.0_R8
          write(*,*) '"ztrgeq" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='ztrgeq'
          call set_targetbarlabel(target_label(itwrk))
        case ('ztrggt')  ! ztrggt
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp = istack(iistack)
          if(itmp.gt.ndepth_nml .or. itmp.le.0) then
            write(*,*) 'No such wall depth index'
            return
          endif
          trg_tmp(1:nbnd,1:itmp) = 0.0_R8
          trg_tmp(1:nbnd,itmp+1:ndepth_nml) = 1.0_R8
          write(*,*) '"ztrggt" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='ztrggt'
          call set_targetbarlabel(target_label(itwrk))
        case('trgm+')      ! target stack addition
          stack_req=(/0,0,0,0,0,2/)
          if(bad_args(stack_req)) return
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              target_stack(ibnd,i,itwrk+1)=
     1         target_stack(ibnd,i,itwrk-1)+
     2         target_stack(ibnd,i,itwrk)
            end do
          end do
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              target_stack(ibnd,i,itwrk-1)=
     1         target_stack(ibnd,i,itwrk+1)
            end do
          end do
          target_label(itwrk-1)='Special'
          itwrk=itwrk-1
          call set_targetbarlabel(target_label(itwrk))
        case('trgm-')      ! target stack subtraction
          stack_req=(/0,0,0,0,0,2/)
          if(bad_args(stack_req)) return   ! 5900
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              target_stack(ibnd,i,itwrk+1)=
     1         target_stack(ibnd,i,itwrk-1)-
     2         target_stack(ibnd,i,itwrk)
            end do
          end do
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              target_stack(ibnd,i,itwrk-1)=target_stack(ibnd,i,itwrk+1)
            end do
          end do
          target_label(itwrk-1)='Special'
          itwrk=itwrk-1
          call set_targetbarlabel(target_label(itwrk))
        case('trgm*')      ! target stack multiplication
          stack_req=(/0,0,0,0,0,2/)
          if(bad_args(stack_req)) return   ! 5900
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              target_stack(ibnd,i,itwrk+1)=
     1         target_stack(ibnd,i,itwrk-1)*
     2         target_stack(ibnd,i,itwrk)
            end do
          end do
          do ibnd=1,nbnd
            do i=1,ndepth_nml
               target_stack(ibnd,i,itwrk-1)=target_stack(ibnd,i,itwrk+1)
            end do
          end do
          target_label(itwrk-1)='Special'
          itwrk=itwrk-1
          call set_targetbarlabel(target_label(itwrk))
        case('trgm/')      ! target stack division
          stack_req=(/0,0,0,0,0,2/)
          if(bad_args(stack_req)) return   ! 5900
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              if(target_stack(ibnd,i,itwrk) .ne. 0.0_R8) then
               target_stack(ibnd,i,itwrk+1)=
     1          target_stack(ibnd,i,itwrk-1)/
     2          target_stack(ibnd,i,itwrk)
              else
                target_stack(ibnd,i,itwrk+1)=0.0_R8
              end if
            end do
          end do
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              target_stack(ibnd,i,itwrk-1)=target_stack(ibnd,i,itwrk+1)
            end do
          end do
          target_label(itwrk-1)='Special'
          itwrk=itwrk-1
          call set_targetbarlabel(target_label(itwrk))
        case default TFVALUE
          cflag=.false.
      END SELECT TFVALUE
c
 6900 continue
      if(cflag) then
        if(b2plot_debug.ge.4)
     1   write(*,*) 'Applied operator ',trim(cval),
     2   ' to target stack'
        if(itwrk.gt.0) then
          values(32)=minval(target_stack(1:nbnd,1:ndepth_nml,itwrk))
          values(33)=maxval(target_stack(1:nbnd,1:ndepth_nml,itwrk))
        endif
      endif
      return

      end subroutine select_target_stack

      subroutine select_constants
      character*7 stack_name

      cflag=.true.
      CNVALUE: select case(cval)
        case('me')      ! electron mass
          stack_name='real'
          call rpush(rstack,irstack,me)
        case('mp')      ! proton mass
          stack_name='real'
          call rpush(rstack,irstack,mp)
        case('ev')      ! electron charge
          stack_name='real'
          call rpush(rstack,irstack,ev)
        case('qe')      ! electron charge
          stack_name='real'
          call rpush(rstack,irstack,qe)
        case('mu0')     ! magnetic permittivity of the vacuum
          stack_name='real'
          call rpush(rstack,irstack,mu0)
        case('eps0')    ! electric permeability of the vacuum
          stack_name='real'
          call rpush(rstack,irstack,eps0)
        case('pi')      ! pi
          stack_name='real'
          call rpush(rstack,irstack,pi)
        case('c')       ! speed of light
          stack_name='real'
          call rpush(rstack,irstack,c)
        case('avogr')   ! Avogadro's number
          stack_name='real'
          call rpush(rstack,irstack,avogr)
        case('kbolt')   ! Boltzmann constant
          stack_name='real'
          call rpush(rstack,irstack,kbolt)
        case('jxi')   ! jxi
          stack_name='integer'
          call ipush(istack,iistack,jxi)
        case('jxa')   ! jxa
          stack_name='integer'
          call ipush(istack,iistack,jxa)
        case('jsep')   ! jsep
          stack_name='integer'
          call ipush(istack,iistack,jsep)
        case('nx')   ! nx 
          stack_name='integer'
          call ipush(istack,iistack,nx)
        case('ny')   ! ny 
          stack_name='integer'
          call ipush(istack,iistack,ny)
        case('ns')   ! ns 
          stack_name='integer'
          call ipush(istack,iistack,ns)
        case('nspecies')   ! nspecies 
          stack_name='integer'
          call ipush(istack,iistack,nspecies)
        case('nxtl')   ! nxtl 
          stack_name='integer'
          if (nnreg(0).eq.8) then
            call ipush(istack,iistack,nxtl)
          else
            write(*,*) cval,
     .       'option available for DN cases only: Nothing done.'
          endif
        case('nxtr')   ! nxtr 
          stack_name='integer'
          if (nnreg(0).eq.8) then
            call ipush(istack,iistack,nxtr)
          else
            write(*,*) cval,
     .       'option available for DN cases only: Nothing done.'
          endif
        case('nbnd')   ! nbnd
          stack_name='integer'
          call ipush(istack,iistack,nbnd)
        case('ndepth')   ! ndepth_nml
          stack_name='integer'
          call ipush(istack,iistack,ndepth_nml)
        case('ismain')   ! ismain
          stack_name='integer'
          call ipush(istack,iistack,ismain)
        case default CNVALUE
          cflag=.false.
      END SELECT CNVALUE

      if(cflag) write(*,*) 'Copied constant ',trim(cval),' to ',
     . trim(stack_name),' stack.'
      return

      end subroutine select_constants

      subroutine select_flvalue

      cflag=.true.
      FLVALUE: select case(cval)
        case('f30')     ! write fort.30
          call write_f30(nx, ny)
        case('f31')     ! write fort.31
          call write_f31(nx, ny)
        case('glob')    ! transform an extra label into a global header
          call set_globalheader
        case('both')    ! both phys and comp
          typeofplot=3  ! both
        case('phys')    ! physical domain only
          typeofplot=1  ! physical domain only
        case('comp')    ! computational domain only
          typeofplot=2  ! computational domain only
        case('mesh')    ! overlay mesh
          call mesh
        case('page')    ! set number of graphs per page
          stack_req=(/2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call b2page(istack(iistack-1),istack(iistack))
          call pop_stack(stack_req)
        case('logf')    ! log f plot
          if (mod(linlog,2).eq.1) linlog=linlog+1
        case('linf')    ! linear f plot
          if (linlog.eq.2.or.linlog.eq.4) linlog=linlog-1
        case('logx')    ! log x plot
          if(linlog.le.2) linlog=linlog+2
        case('linx')    ! linear x plot
          if(linlog.gt.2) linlog=linlog-2
        case('clvs')    ! Customise Contour Levels
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          stack_req=(/1,istack(iistack),0,0,0,0/)
          if(bad_args(stack_req)) return
          call clvs(istack(iistack),rstack,irstack)
          do iclv=1,istack(iistack)
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('clvd')    ! Customise Contour Level Delta
          stack_req=(/0,3,0,0,0,0/)
          if(bad_args(stack_req)) return
          call clvd(rstack(irstack-2),rstack(irstack-1),rstack(irstack))
          do iclv=1,3
            dummy=rpop(rstack,irstack)
          end do
        case('clvn')    ! Customise Contour Level Range/Num
          stack_req=(/1,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          call clvn(rstack(irstack-1),rstack(irstack),istack(iistack))
          do iclv=1,2
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('clvp')    ! Customise Contour Level by Proportion
          stack_req=(/2,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          ix=ipop(istack,iistack)
          iy=ipop(istack,iistack)
          call clvp(ix,iy)
          write(stmp,'(a,2x,a1,2i4,a2)') 
     -     graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk))),'{',ix,iy,' }'
          graphlabel(iwrk)=stmp
        case('debu','debug')  ! set debugging level
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          b2plot_debug=ipop(istack,iistack)
          write(*,*) 'Debugging set to ',b2plot_debug
        case('eclv')    ! Do not fill Contour Levels
          call nofill(.true.)
        case('fclv')    ! Do fill Contour Levels
          call nofill(.false.)
        case('font')    ! Set font
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call font(istack(iistack))
          call pop_stack(stack_req)
        case('lacl')    ! Lbl Orient. along Contour Levels
          call lblalongcl(.true.)
        case('lhor')    ! Label Orientation horizontal
          call lblalongcl(.false.)
        case('lbon')    ! Labeling on / off
          call labelon(.true.)
        case('lbof')    ! Labeling on / off
          call labelon(.false.)
        case('lwid')    ! Line Width Factor
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          call linewidth(rstack(irstack))
          lwidth=rpop(rstack,irstack)
        case('zwid')    ! Axis Line Width Factor
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          call axislinewidth(rstack(irstack))
          zwidth=rpop(rstack,irstack)
        case('lfsz')    ! Line Label Font Size
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_line_font_size(istack(iistack))
          lfsz=istack(iistack)
          call pop_stack(stack_req)
        case('grid')    ! Grid on
          call gridon(.true.)
        case('clvl')    ! Contour Level Lines off
          call clvl(.false.)
        case('spln')    ! switch on interpolation
          call ipolate(.true.)
        case('arow')    ! change default arrow length
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          arrowl=rpop(rstack,irstack)
        case('arww')    ! change default arrow width
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          arroww=rpop(rstack,irstack)
        case('acol')    ! colour index for arrows
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          acol=istack(iistack)
          call pop_stack(stack_req)
        case('sep')     ! plot separatrix
          call separatrixon(.true.)
        case('vesl')    ! plot vessel etc.
          call vesselon(.true.)
        case('mmid')    ! map to midplane
          call set_midplane_mapping(.true.)
        case('nocp')    ! do not plot contour plot
          call contourplot(.false.)
        case('sepw')    ! change separatrix line width
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          sepwidth=rpop(rstack,irstack)
        case('ewid')    ! change experimental data line width
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          ewidth=rpop(rstack,irstack)
        case('owid')    ! change outline line width
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          owidth=rpop(rstack,irstack)
        case('vwid')    ! change vessel line width
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          vwidth=rpop(rstack,irstack)
        case('sq')      ! square aspect ratio
          call set_square
        case('xpol')      ! use poloidal distance for f.x
         xparallel=.false.
        case('xpar')      ! use parallel distance for f.x
         xparallel=.true.
        case('a4p')     ! A4 portrait
          call set_a4p
        case('a4l')     ! A4 landscape
          call set_a4l
        case('same')    ! next graph on same axes
          call set_same_page
        case('filo')    ! fill option
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_fill_option(istack(iistack))
          call pop_stack(stack_req)
        case('i+')      ! i2 + i1
          stack_req=(/2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
	  it2=ipop(istack,iistack)
	  it1=ipop(istack,iistack)
          call ipush(istack,iistack,it1+it2)
        case('i-')      ! i2 - i1
          stack_req=(/2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
	  it2=ipop(istack,iistack)
	  it1=ipop(istack,iistack)
          call ipush(istack,iistack,it1-it2)
        case('i*')      ! i2 * i1
          stack_req=(/2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
	  it2=ipop(istack,iistack)
	  it1=ipop(istack,iistack)
          call ipush(istack,iistack,it1*it2)
        case('i/')      ! i2 / i1
          stack_req=(/2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          it2=ipop(istack,iistack)
          it1=ipop(istack,iistack)
          if (it2.eq.0) then
            write(*,*) 'Division by zero : Operation not permitted'
            return
          endif
          call ipush(istack,iistack,it1/it2)
        case('i**')     ! i2 ** i1
          stack_req=(/2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          it2=ipop(istack,iistack)
          it1=ipop(istack,iistack)
          if (it1.eq.0 .and. it2.eq.0) then
            write(*,*) 'Zero to the zeroth power : Undefined'
            return
          endif
          call ipush(istack,iistack,it1**it2)
        case('itor')    ! integer --> real stack
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          dummy=real(ipop(istack,iistack))
          call rpush(rstack,irstack,dummy)
        case('itostr')  ! integer --> string stack
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          write(hlp_frm,'(a,i2,a)') 
     .     '(i',len_of_digits(istack(iistack)),')'
          write(stmp,hlp_frm) ipop(istack,iistack)
          call spush(sstack,isstack,stmp)
        case('itom')    ! integer --> matrix stack
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(iwrk .eq. wklng) call dropstackbottom()
          iwrk=iwrk+1
          ndim(1,iwrk)=nx
          ndim(2,iwrk)=ny
          ndim(3,iwrk)=1
          work(-1:nx,-1:ny,0,iwrk)=real(istack(iistack))
          values(9)=real(istack(iistack))
          values(10)=real(istack(iistack))
          textwk(0,iwrk)=textbl(0)
          write(hlp_frm,'(a,i2,a)') 
     .     '(i',len_of_digits(istack(iistack)),')'
          write(graphlabel(iwrk),hlp_frm) ipop(istack,iistack)
        case('itow')    ! integer --> wall stack
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(iwwrk .eq. wklng) call dropwallstackbottom()
          iwwrk=iwwrk+1
          nwdim(1,iwwrk)=nbnd
          nwdim(2,iwwrk)=1
          wall_stack(1:nbnd,0,iwwrk)=real(istack(iistack))
          values(14)=real(istack(iistack))
          values(15)=real(istack(iistack))
          textwwk(0,iwwrk)=textbl(0)
          write(hlp_frm,'(a,i2,a)') 
     .     '(i',len_of_digits(istack(iistack)),')'
          write(wall_label(iwwrk),hlp_frm) ipop(istack,iistack)
          call set_wallbarlabel(wall_label(iwwrk))
        case('itot')    ! integer --> target stack
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(itwrk .eq. wklng) call droptargetstackbottom()
          itwrk=itwrk+1
          target_stack(1:nbnd,1:ndepth_nml,itwrk)=real(istack(iistack))
          values(32)=real(istack(iistack))
          values(33)=real(istack(iistack))
          write(hlp_frm,'(a,i2,a)') 
     .     '(i',len_of_digits(istack(iistack)),')'
          write(target_label(itwrk),hlp_frm) ipop(istack,iistack)
          call set_targetbarlabel(target_label(itwrk))
        case('iprt')    ! output the integer stack
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is = iistack, 1, -1
            write(*,*) istack(is)
            if(writing_to_file) then
              write(nwtf,*) istack(is)                   !srv 21.09.11
            endif
          enddo
        case('r+')      ! r2 + r1
          stack_req=(/0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
	  t2=rpop(rstack,irstack)
	  t1=rpop(rstack,irstack)
          call rpush(rstack,irstack,t1+t2)
        case('r-')      ! r2 - r1
          stack_req=(/0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
	  t2=rpop(rstack,irstack)
	  t1=rpop(rstack,irstack)
          call rpush(rstack,irstack,t1-t2)
        case('r*')      ! r2 * r1
          stack_req=(/0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
	  t2=rpop(rstack,irstack)
	  t1=rpop(rstack,irstack)
          call rpush(rstack,irstack,t1*t2)
        case('r/')      ! r2 / r1
          stack_req=(/0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          t2=rpop(rstack,irstack)
          t1=rpop(rstack,irstack)
          if (t2.eq.0.0_R8) then
            write(*,*) 'Division by zero : Operation not permitted'
            return
          endif
          call rpush(rstack,irstack,t1/t2)
        case('r**')     ! r2 ** r1
          stack_req=(/0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          t2=rpop(rstack,irstack)
          t1=rpop(rstack,irstack)
          if (t1.eq.0.0_R8 .and. t2.eq.0.0_R8) then
            write(*,*) 'Zero to the zeroth power : Undefined'
            return
          endif
          call rpush(rstack,irstack,t1**t2)
        case('rtostr')  ! real --> string stack
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          write(stmp,'(1pe15.8)') rpop(rstack,irstack)
          call spush(sstack,isstack,stmp)
        case('rprt')    ! output the real stack
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is = irstack, 1, -1
            write(*,'(1pe15.8)') rstack(is)
            if(writing_to_file) then
              write(nwtf,'(1pe15.8)') rstack(is)                   !srv 21.09.11
            endif
          enddo
        case('new')     ! new page
          call set_new_page
        case('lbsz')    ! label size
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_lbsz(rpop(rstack,irstack))
        case('ndec')    ! number of decades
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_ndec(istack(iistack))
          call pop_stack(stack_req)
        case('aspc')    ! preserve aspect ratio flag
          call swap_aspect
        case('col')     ! colour
#ifndef PLPLOT
          call sfseti('TY',0)
#endif
          color_set=.true.
        case('bw')      ! black and white
#ifndef PLPLOT
          call sfseti('TY',-4)
#endif
          color_set=.false.
        case('logt')    ! log time for trajectories
          tlinlog=0
        case('lint')    ! linear time for trajectories
          tlinlog=1
        case('zero')    ! plot origin
          call set_zero
        case('outl')    ! outline of plasma grid
          call outlineon(.true.)
        case('stdl')    ! "standard" label
          call customise_label
     .     ('${HOST} ${PWD} `date` T=${TIME}',stmp)
          call clear_extralabel(stmp)
        case('mvsk')    ! movie skip
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          movie_skip=istack(iistack)
          call pop_stack(stack_req)
        case('ave')     ! number of things to average over
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          average=istack(iistack)
          call pop_stack(stack_req)
        case('t_#p')    ! number of trajectories
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call t_setnplot(istack(iistack))
          call pop_stack(stack_req)
        case('reso')    ! reset output
          call resetoutput
        case('swid')    ! part of tcl/tk interface
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call setwinid(istack(iistack))
          call pop_stack(stack_req)
        case('xyrz')    ! part of tcl/tk interface
          call getxyrz
        case('t_r+')    ! ratio true
          ratio=.true.
        case('t_r-')    ! ratio false
          ratio=.false.
        case('mcol')    ! colour index of mesh
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          mcol=istack(iistack)
          call pop_stack(stack_req)
        case('scol')    ! colour index of separatrix
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          scol=istack(iistack)
          call pop_stack(stack_req)
        case('vcol')    ! colour index of vessel
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          vcol=istack(iistack)
          call pop_stack(stack_req)
        case('ocol')    ! colour index of grid outline
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          ocol=istack(iistack)
          call pop_stack(stack_req)
        case('zcol')    ! colour index of zero lines
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_axes_color(istack(iistack))
          zcol=istack(iistack)
          call pop_stack(stack_req)
        case('hcol')    ! colour index of text header
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_header_color(istack(iistack))
          hcol=istack(iistack)
          call pop_stack(stack_req)
        case('gcol')    ! colour index of graph title
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_title_color(istack(iistack))
          gcol=istack(iistack)
          call pop_stack(stack_req)
        case('bcol')    ! colour index of bar label title
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_barlabel_color(istack(iistack))
          bcol=istack(iistack)
          call pop_stack(stack_req)
        case('xcol')    ! colour index of extralabel
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_extralabel_color(istack(iistack))
          xcol=istack(iistack)
          call pop_stack(stack_req)
        case('dcol')    ! colour index for matrix data line plots
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          dcol=istack(iistack)
          call pop_stack(stack_req)
        case('wcol')    ! colour index for wall data line plots
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wcol=istack(iistack)
          call pop_stack(stack_req)
        case('tcol')    ! colour index for target data line plots
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          tcol=istack(iistack)
          call pop_stack(stack_req)
        case('rigz')    ! toggle set\_resignore\_regions\_to\_zero
          set_resignore_regions_to_zero=
     1     .not.set_resignore_regions_to_zero  
        case('noax')    ! no axes
          call set_noax
        case('part')    ! partition colour map
          stack_req=(/0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_partition_colors(rpop(rstack,irstack))
        case('lbsa')    ! toggle label lines
          call toggle_label_lines
        case('ytra')    ! transform the "y" coordinate for comparison to experiment
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          y_transform_l=istack(iistack)
          call pop_stack(stack_req)
        case('ecol')    ! colour index of experimental data
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          ecol=istack(iistack)
          call pop_stack(stack_req)
        case('nspe','nspec')    ! new spectrometer
          stack_req=(/3,0,0,0,0,0/)
          if(bad_args(stack_req)) return
c     Parameters on istack: <species number> <wavelength in A>
c     <delta wavelength in A> nspe
          call CalculateEmissivity(istack(iistack-2),
     &     istack(iistack-1),
     &     istack(iistack))
          call pop_stack(stack_req)
        case('spec')    ! old spectrometer
          call spect
        case('coro')    ! print coronal rates
          stack_req=(/1,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          call callgetsa(ns,rstack(irstack),istack(iistack))
          call pop_stack(stack_req)
        case('cora')    ! coronal radiation rate
          call setcorora(ismain)
        case('coni')    ! coronal densities
          call setcoroni(ismain)
        case('shco')    ! part of tcl/tk
          continue
        case('user')    ! call a user function "b2paux"
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call b2paux(istack(iistack))
          call pop_stack(stack_req)
        case('writ','write')    ! set writing to file
          if(.not.writing_to_file) then
            writing_to_file=.true.
            open(nwtf,file='b2plot.write')                           !srv 21.09.11
          end if
        case('omp')     ! set jxa (OMP)
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.nx) then
            jxa = ipop(istack,iistack)
            write(*,*) 'jxa changed to ', jxa
          else
            write(*,*) 
     .       'Faulty ix = ',istack(iistack),' chosen for jxa!'
            call pop_stack(stack_req)
          endif
        case('imp')     ! set jxi (IMP)
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.nx) then
            jxi = ipop(istack,iistack)
            write(*,*) 'jxi changed to ', jxi
          else
            write(*,*) 
     .       'Faulty ix = ',istack(iistack),' chosen for jxi!'
            call pop_stack(stack_req)
          endif
        case('yref')     ! set yref
          stack_req=(/1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(-2.le.istack(iistack).and.istack(iistack).le.ny) then
            iyref = ipop(istack,iistack)
            if (iyref.gt.-2) then
              write(*,*) 'iyref changed to ', iyref
            else
              write(*,*) 'iyref option turned off'
            endif
          else
            write(*,*) 
     .       'Faulty iy = ',istack(iistack),' chosen for iyref!'
            call pop_stack(stack_req)
          endif
        case default FLVALUE
          cflag=.false.
      END SELECT FLVALUE

      if (cflag.and.b2plot_debug.ge.1) write(*,*) trim(cval),' : Done.'
      return

      end subroutine select_flvalue

      LOGICAL FUNCTION bad_args(stack_req)
      integer stack_req(6)
      bad_args=.false.
      if(stack_req(1).gt.iistack) then
        write(*,*) 'Insufficient data on integer stack'
        bad_args= .true.
      endif
      if(stack_req(2).gt.irstack) then
        write(*,*) 'Insufficient data on real stack'
        bad_args= .true.
      endif
      if(stack_req(3).gt.isstack) then
        write(*,*) 'Insufficient data on string stack'
        bad_args= .true.
      endif
      if(stack_req(4).gt.iwrk) then
        write(*,*) 'Insufficient data on matrix stack'
        bad_args= .true.
      endif
      if(stack_req(5).gt.iwwrk) then
        write(*,*) 'Insufficient data on wall stack'
        bad_args= .true.
      endif
      if(stack_req(6).gt.itwrk) then
        write(*,*) 'Insufficient data on target stack'
        bad_args= .true.
      endif
      END FUNCTION bad_args

      subroutine pop_stack(stack_req)
      integer stack_req(6),is,itmp
      do is=1,stack_req(1)
        itmp=ipop(istack,iistack)
      enddo
      end subroutine pop_stack

      end subroutine b2plot

      subroutine write_crazy(lu,ndata,kind,is,s_ref,q)
      use b2mod_types
      implicit none
      integer lu,ndata,kind,is,i
      real(kind=R8) s_ref(ndata),q(ndata)

      write(lu,'(1x,a,i2,a,i3,a,i4,a)') 
     .    'ndata(1,',kind,',',is,')=',ndata,','
      do i=1,ndata
        write(lu,'(1x,1p,2(a,i4,a,i2,a,i3,a,1e15.8,a))') 
     .    'tdata(1,',i,',',kind,',',is,')=',s_ref(i),',',
     1   ' tdata(2,',i,',',kind,',',is,')=',q(i),','
      enddo
      end subroutine write_crazy

!#ifdef USE_EIRENE_NOT_YET_READY
      subroutine get_amjuel_data(filnam,h123,reac,crc,
     1 nxd,nyd,nx,ny,ne,te,mach)
      use b2mod_types
      use b2mod_constants
      IMPLICIT NONE
#     include <eir_reac.F>
      integer nxd,nyd,nx,ny
      real (kind=R8) ::
     1 ne(0:nxd+1,0:nyd+1),te(0:nxd+1,0:nyd+1),
     2 mach(0:nxd+1,0:nyd+1),D11(0:8,0:8)
      CHARACTER FILNAM*8,H123*4,REAC*9,CRC*3
      integer ir,i,j,ix,iy
      real (kind=R8) ::
     1 rmin,rmax,tef,def,dnnn,tei,dej
      IR=1
      CALL SLREAC_b2plot(IR,FILNAM,H123,REAC,CRC)
      DO J=1,9
        DO I=1,9
          if(h123.eq.'H.-1') then
            D11(J-1,I-1)=CREAC(J,0,IR)
          elseif(
     1     h123.eq.'H.1'.or.
     1     h123.eq.'H.2'.or.
     1     h123.eq.'H.5'.or.
     1     h123.eq.'H.8'.or.
     1     h123.eq.'H.11') then
            D11(J-1,I-1)=CREAC(J,1,IR)
          elseif(
     1     h123.eq.'H.3'.or.
     1     h123.eq.'H.4'.or.
     1     h123.eq.'H.6'.or.
     1     h123.eq.'H.7'.or.
     1     h123.eq.'H.9'.or.
     1     h123.eq.'H.10'.or.
     1     h123.eq.'H.12') then
            D11(J-1,I-1)=CREAC(J,I,IR)
          else
            write(*,*) 'Don''t yet know about ',h123
            D11(J-1,I-1)=CREAC(J,I,IR)
          endif
        enddo
      enddo

      rmin=1.0e30
      rmax=1.0e-30
      do iy=0,ny+1
        do ix=0,nx+1
          TEF=LOG(TE(ix,iy)/qe)
          DEF=LOG(ne(ix,iy)*1.0d-8*1.0d-6)
          dnnn=0.0
          DO I=0,8
            TEI=TEF**I
            if(
     1       h123.eq.'H.1'.or.
     1       h123.eq.'H.2'.or.
     1       h123.eq.'H.5'.or.
     1       h123.eq.'H.8'.or.
     1       h123.eq.'H.11') then
              DNNN=DNNN+D11(I,0)*TEI
            elseif(
     1       h123.eq.'H.3'.or.
     1       h123.eq.'H.4'.or.
     1       h123.eq.'H.6'.or.
     1       h123.eq.'H.7'.or.
     1       h123.eq.'H.9'.or.
     1       h123.eq.'H.10'.or.
     1       h123.eq.'H.12') then
              DO J=0,8
                DEJ=DEF**J
                DNNN=DNNN+D11(I,J)*TEI*DEJ
              enddo
            else
c             write(*,*) 'Don''t yet know about ',h123
            endif
          enddo
          if(dnnn.gt.log(1.0e30)) then
            dnnn=1.0e30
          elseif(dnnn.lt.log(1.0e-30)) then
            dnnn=1.0e-30
          else
            DNNN=EXP(DNNN)
          endif
          rmin=min(rmin,dnnn)
          rmax=max(rmax,dnnn)
          mach(ix,iy)=dnnn
        enddo
      enddo
      write(*,*) 'rmin,rmax= ',rmin,rmax
      return
      end subroutine get_amjuel_data

      SUBROUTINE SLREAC_b2plot (IR,FILNAM,H123,REAC,CRC)
      use b2mod_types
c
C  input
C    filnam: read a&m data from file filnam, e.g. AMJUEL
c    ir    : store data on eirene array creac(...,...,ir)
c    H123  : identifyer for reaction in filnam, e.g. H2
c    REAC  : number of reaction in filnam, e.g. 2.2.5
c    crc   : type of data process e.g. EI, CX, OT, etc
C  output
c    iswr  : eirene flag for type of process  (1,2,...7)
c    creac : eirene storage array for a&m data CREAC(9,0:9,IR)
c    modcol: see below
c
      IMPLICIT NONE
C  READ A&M DATA FROM THE FILES INTO EIRENE ARRAY CREAC
C
C  INPUT (PARAMETER LIST):
C    IR: REACTION NUMBER IN CREAC
C    FILNAM: HYDHEL, METHANE, AMJUEL OR CONST.
C
C  OUTPUT (IN COMMON COMXS):
C    READ DATA FROM "FILNAM" INTO ARRAY "CREAC"
C    DEFINE PARAMETER MODCLF(IR) (4 DIGITS NMLK)
C    FIRST DECIMAL  K           =1  CROSS SECTION AVAILABLE
C                                   (ON CREAC(..,0,IR))
C                   K           =0  ELSE
C    SECOND DECIMAL L           =1  <SIGMA V> FOR ONE
C                                   PARAMETER E (E.G.
C                                   PROJECTILE ENERGY OR ELECTRON
C                                   DENSITY) AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  <SIGMA V> FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  <SIGMA V> FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   L           =0  ELSE
C    THIRD  DECIMAL M               DATA FOR MOMENTUM EXCHANGE
C                                   TO BE WRITTEN
C    FOURTH DECIMAL N           =1  DELTA E FOR ONE (E.G. E=0)
C                                   PROJECTILE ENERGY AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  DELTA E FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  DELTA E FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   N           =0  ELSE
C
#      include <eir_reac.F>
      real (kind=R8) ::
     1 CONST
      integer IR,I0,ISW,J,IND,I,K,IH,ireac
      CHARACTER ZEILE*80
      CHARACTER FILNAM*8,H123*4,REAC*9,CRC*3
      CHARACTER*2 CHR
      CHARACTER*3 CHRL,CHRR
      LOGICAL LCONST,LGEMIN,LGEMAX
      logical open_file
      external open_file
C
      if(ir.gt.nreac) stop 'SLREAC_b2plot: IR > NREAC'
      LGEMIN=.FALSE.
      LGEMAX=.FALSE.
      ISWR(IR)=0
      CONST=0.
      CHR='l0'
      I0=0
C
      IF (INDEX(CRC,'EI').NE.0) ISWR(IR)=1
      IF (INDEX(CRC,'DS').NE.0) ISWR(IR)=2
      IF (INDEX(CRC,'CX').NE.0) ISWR(IR)=3
      IF (INDEX(CRC,'II').NE.0) ISWR(IR)=4
      IF (INDEX(CRC,'EL').NE.0) ISWR(IR)=5
      IF (INDEX(CRC,'RC').NE.0) ISWR(IR)=6
      IF (INDEX(CRC,'OT').NE.0) ISWR(IR)=7
C
      IF (INDEX(FILNAM,'AMJUEL').NE.0) THEN
        LCONST=.FALSE.
        if(.not.open_file(29,'AMJUEL',
     1   'data.local|data')) stop 'No AMJUEL found'
      ELSEIF (INDEX(FILNAM,'METHAN').NE.0) THEN
        LCONST=.FALSE.
        if(.not.open_file(29,'METHANE',
     1   'data.local|data')) stop 'No METHANE found'
      ELSEIF (INDEX(FILNAM,'HYDHEL').NE.0) THEN
        LCONST=.FALSE.
        if(.not.open_file(29,'HYDHEL',
     1   'data.local|data')) stop 'No HYDHEL found'
      ELSE
        LCONST=.TRUE.
        READ (REAC(1:9),'(E9.2)') CONST
        IF (CONST.LE.0.) GOTO 991
      ENDIF
C
      IF (H123(4:4).EQ.' ') THEN
        READ (H123(3:3),'(I1)') ISW
      ELSE
        READ (H123(3:4),'(I2)') ISW
      ENDIF
C
      IREAC=INDEX(REAC,' ')-1
      IF (IREAC.LT.0) IREAC=LEN(REAC)

C
      IF (ISW.EQ.1) THEN
        CHR='a0'
        CHRL='al0'
        CHRR='ar0'
        I0=0
        MODCLF(IR)=MODCLF(IR)+1
      ELSEIF (ISW.EQ.2) THEN
        CHR='b0'
        CHRL='bl0'
        CHRR='br0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+10
      ELSEIF (ISW.EQ.3) THEN
        MODCLF(IR)=MODCLF(IR)+20
      ELSEIF (ISW.EQ.4) THEN
        MODCLF(IR)=MODCLF(IR)+30
      ELSEIF (ISW.EQ.5) THEN
        CHR='e0'
        CHRL='el0'
        CHRR='er0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+100
      ELSEIF (ISW.EQ.6) THEN
        MODCLF(IR)=MODCLF(IR)+200
      ELSEIF (ISW.EQ.7) THEN
        MODCLF(IR)=MODCLF(IR)+300
      ELSEIF (ISW.EQ.8) THEN
        CHR='h0'
        CHRL='hl0'
        CHRR='hr0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+1000
      ELSEIF (ISW.EQ.9) THEN
        MODCLF(IR)=MODCLF(IR)+2000
      ELSEIF (ISW.EQ.10) THEN
        MODCLF(IR)=MODCLF(IR)+3000
      ELSEIF (ISW.EQ.11) THEN
        CHR='k0'
        CHRL='kl0'
        CHRR='kr0'
        I0=1
      ELSEIF (ISW.EQ.12) THEN
      ENDIF
C
      IF (.NOT.LCONST) THEN
100     READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'##BEGIN DATA HERE##').EQ.0) GOTO 100

1       READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,H123).EQ.0) GOTO 1
C
2       READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'H.').NE.0) GOTO 990
        IF (INDEX(ZEILE,'Reaction ').EQ.0.or.
     .      INDEX(ZEILE,REAC(1:ireac)).EQ.0) GOTO 2
      ENDIF
C
C  SINGLE PARAM. FIT, ISW=1,2,5,8,11
      IF (ISW.EQ.1.OR.ISW.EQ.2.OR.ISW.EQ.5.OR.ISW.EQ.8.OR.ISW.EQ.11)
     .THEN
        IF (.NOT.LCONST) THEN
3         READ (29,'(A80)',END=990) ZEILE
          IF (INDEX(ZEILE,CHR).EQ.0) GOTO 3
          DO 9 J=0,2
            IND=0
            DO 4 I=1,3
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
              READ (ZEILE((IND+2):80),'(E20.12)') CREAC(J*3+I,I0,IR)
4           CONTINUE
            READ (29,'(A80)',END=990) ZEILE
9         CONTINUE
          IF (INDEX(ZEILE,CHRL).NE.0) THEN
            IND=0
            DO 5 I=1,3
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
              READ (ZEILE((IND+3):80),'(E20.12)') fparm(IR,i)
5           CONTINUE
            LGEMIN=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
          IF (INDEX(ZEILE,CHRR).NE.0) THEN
            IND=0
            DO 7 I=4,6
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
              READ (ZEILE((IND+3):80),'(E20.12)') fparm(IR,i)
7           CONTINUE
            LGEMAX=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
c
          rcmn(ir)=-20.
          rcmx(ir)= 20.
          ifexmn(ir)=0
          ifexmx(ir)=0
c
          if (lgemin) then
            IND=INDEX(ZEILE,'=')
            READ (ZEILE((IND+2):80),'(E12.5)') rcmn(IR)
            rcmn(ir)=log(rcmn(ir))
            ifexmn(ir)=5
            READ (29,'(A80)',END=990) ZEILE
          endif
          if (lgemax) then
            IND=INDEX(ZEILE,'=')
            READ (ZEILE((IND+2):80),'(E12.5)') rcmx(IR)
            rcmx(ir)=log(rcmx(ir))
            ifexmx(ir)=5
            READ (29,'(A80)',END=990) ZEILE
          endif
C
        ELSE
          CREAC(1,I0,IR)=LOG(CONST)
          DO 6 J=2,9
            CREAC(J,I0,IR)=0.
6         CONTINUE
        ENDIF
C
      ELSE
C  TWO PARAM. FIT, ISW=3,4,6,7,9,10,12

        DO 11 J=0,2
16        READ (29,'(A80)',END=990) ZEILE
          IF (INDEX(ZEILE,'Index').EQ.0) GOTO 16
          READ (29,'(1X)')
          DO 17 I=1,9
            READ (29,*) IH,(CREAC(I,K,IR),K=J*3+1,J*3+3)
cdr         write (6,*) IH,(CREAC(I,K,IR),K=J*3+1,J*3+3)
17        CONTINUE
11      CONTINUE
C
      ENDIF
C
      CLOSE (UNIT=29)
C
      RETURN
C
990   WRITE (6,*) ' NO DATA FOUND FOR REACTION ',H123,' ',REAC,
     .            ' IN DATA SET ',FILNAM
      WRITE (6,*) ' IR,MODCLF(IR) ',IR,MODCLF(IR)
      CLOSE (UNIT=29)
      RETURN
991   WRITE (6,*) ' INVALID CONSTANT IN SLREAC. CONST= ',CONST
      WRITE (6,*) ' CHECK "REACTION CARDS" FOR REACTION NO. ',IR
      STOP
      END SUBROUTINE SLREAC_b2plot
!#endif

      subroutine old_chang_hinton(nx,ny,ns,ismain,jsep,jxa,lambda,
     1 safetyfactor,cr,cz,crx,bb,ne,te,ti,na,transport)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nx,ny,ns,ismain,jsep,jxa
      real (kind=R8) :: lambda,
     1 cr(-1:nx,-1:ny), cz(-1:nx,-1:ny), crx(-1:nx,-1:ny,0:3),
     2 bb(-1:nx,-1:ny,0:3),
     2 ne(-1:nx,-1:ny), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     3 na(-1:nx,-1:ny,0:ns-1), safetyfactor(-1:ny),
     2 transport(-1:nx,-1:ny,0:ns-1)

      integer ix,iy,is,irmin,irmax
      real (kind=R8) :: crmin,crmax,czmin,czmax,
     1 cr2max,cr2min,cr1max,cr1min,
     2 bave,epsilon,lnlam,vthi,omegaci,
     1 rhoi,alpha,deltap,f1,f2,taui,nui,nustari,mustari,chi

      transport=0
      do iy=-1,jsep
        write(*,*) 'iy: ',iy
        ix=jxa
        if(mod(region(ix,iy,0),4).ne.1) then
          write(*,*) 'We have a problem: not in core!'
          stop
        endif
        irmin=jxa
        irmax=jxa
        crmin=cr(ix,iy)
        crmax=crmin
        czmin=cz(ix,iy)
        czmax=czmin
        ix=rightix(ix,iy)
        do while(ix.ne.jxa.and.ix.ge.-1.and.ix.le.nx)
          if(cr(ix,iy).lt.crmin) irmin=ix
          if(cr(ix,iy).gt.crmax) irmax=ix
          crmin=min(crmin,cr(ix,iy))
          crmax=max(crmax,cr(ix,iy))
          czmin=min(czmin,cz(ix,iy))
          czmax=max(czmax,cz(ix,iy))
          ix=rightix(ix,iy)
        enddo
        if(ix.ne.jxa) then
          write(*,*) 'We have a problem: core not continuous!'
          stop
        endif
        write(*,*) 'crmin,crmax,czmin,czmax',
     1   crmin,crmax,czmin,czmax
        bave=cr(ix,iy)*2.0_R8/(crmax+crmin)*bb(ix,iy,3)
        write(*,*) 'ne,te,ti,b,q: ',
     1   ne(ix,iy),te(ix,iy)/qe,ti(ix,iy)/qe,bave,safetyfactor(iy)
        epsilon=(crmax-crmin)/(crmax+crmin)
        write(*,*) 'epsilon: ',epsilon
cxpb    Coulomb logarithm formula from Wesson for ion-ion collisions
        if (lambda.lt.0.0_R8) then
          lnlam=max(-lambda,17.3_R8-
     1     log(ne(ix,iy)/1.0e20_R8)/2.0_R8+
     2     1.5_R8*log(ti(ix,iy)/ev/1000.0_R8))
        else
          lnlam = lambda
        endif
        write(*,*) 'lnlam: ',lnlam
        cr2max=0.5_R8*(crx(irmax,iy,1)+crx(irmax,iy,3))
        cr2min=0.5_R8*(crx(irmin,iy,1)+crx(irmin,iy,3))
        cr1max=0.5_R8*(crx(irmax,iy,0)+crx(irmax,iy,2))
        cr1min=0.5_R8*(crx(irmin,iy,0)+crx(irmin,iy,2))
        write(*,*) 'cr2max,cr2min,cr1max,cr1min: ',
     1   cr2max,cr2min,cr1max,cr1min
        deltap=((cr2max+cr2min)-(cr1max+cr1min))/
     1         ((cr2max-cr2min)-(cr1max-cr1min))
        write(*,*) 'deltap: ',deltap
        f1=(1_R8+1.5_R8*(epsilon**2+epsilon*deltap)+
     1         0.375_R8*epsilon**3*deltap)/(1_R8+0.5_R8*epsilon*deltap)
        write(*,*) 'f1: ',f1
        f2=sqrt(1_R8-epsilon**2)*(1.0_R8+0.5_R8*epsilon*deltap)/
     1         (1_R8+deltap/epsilon*(sqrt(1_R8-epsilon**2)-1.0_R8))
        write(*,*) 'f2: ',f2
        do is=0,ns-1
          if(.not.is_neutral(is)) then
            write(*,*) 'za,am,na: ',rza(ix,iy,is),am(is),na(ix,iy,is)
            vthi=sqrt(ti(ix,iy)/mp/am(is))
            write(*,*) 'vthi: ',vthi
            omegaci=rza(ix,iy,is)*qe*bave/mp/am(is)
            write(*,*) 'omegaci: ',omegaci
            rhoi=sqrt(2.0_R8)*vthi/omegaci
            write(*,*) 'rhoi: ',rhoi
            alpha=na(ix,iy,is)*rz2(ix,iy,is)/
     1            na(ix,iy,ismain)/rz2(ix,iy,ismain)  !!! think about this !!!
            write(*,*) 'alpha: ',alpha
            taui=12.0_R8*pi**1.5_R8*eps0**2*
     1             (am(is)*mp)**0.5_R8*ti(ix,iy)**1.5_R8/
     2             (na(ix,iy,is)*rz2(ix,iy,is)**2*qe**4*lnlam)
            write(*,*) 'taui: ',taui
            nui=1.0_R8/taui
            nustari=nui*(crmax+crmin)/2.0_R8*safetyfactor(iy)/
     1             epsilon**(1.5_R8)/vthi
            write(*,*) 'nustari: ',nustari
            mustari=nustari*(1.0_R8+1.54_R8*alpha)
            write(*,*) 'mustari: ',mustari
            chi=epsilon**(-1.5_R8)*
     1             safetyfactor(iy)**2*rhoi**2/taui*
     2             ((0.66_R8*(1.0_R8+1.54_R8*alpha)+
     3              (1.88_R8*epsilon**0.5_R8-1.54_R8*epsilon)*
     4              (1_R8+3.75_R8*alpha))/
     5              (1.0_R8+1.03_R8*mustari**0.5_R8+0.31*mustari)*f1
     5              +0.59_R8*mustari*epsilon/
     6               (1_R8+0.74_R8*mustari*epsilon**1.5_R8)*
     7               (1_R8+
     8                 1.33_R8*alpha*(1.0_R8+0.60_R8*alpha)/
     9                 (1.0_R8+1.79_R8*alpha))*(f1-f2))
            write(*,*) 'chi: ',chi
            transport(ix,iy,is)=chi
            ix=rightix(ix,iy)
            do while(ix.ne.jxa.and.ix.ge.-1.and.ix.le.nx)
              transport(ix,iy,is)=chi
              ix=rightix(ix,iy)
            enddo
          endif
        enddo
      enddo
      return
      end subroutine old_chang_hinton
