c The following can be enabled (by removing the ! before the @) to
c speed up compilation of b2plot on IBMs where it can take a long time.
c Unfortunately the Juelich compiler complains with error message
c    (E) The STRICT/NOSTRICT option may only be specified for OPT(3).  STRICT is assumed.
c so this optimisation is disabled by default.
c 2002-09-11: logic revised
      subroutine b2plot(nxd,nyd,nsd_local)
#ifdef NO_OPT
#ifdef CRAY
!DIR$ OPTIMIZE (-O0 -hfp0)
#else
!DIR$ NOOPTIMIZE
#endif
#endif
      use b2mod_types
      use b2mod_geo
      use b2mod_work
      use b2mod_time
      use b2mod_diag
      use b2mod_rates
      use b2mod_b2plot
      use b2mod_movies
      use b2mod_plasma
      use b2mod_sources
      use b2mod_sputter
      use b2mod_indirect
      use b2mod_elements
      use b2mod_constants
      use b2mod_residuals
      use b2mod_transport
      use b2mod_geo_corner
      use b2mod_b2plot_xy
      use b2mod_b2plot_ank
      use b2mod_b2plot_help
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_b2plot_wall_loading
      use b2mod_user_namelist
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_anomalous_transport
      use b2mod_neoclassical
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_ranges
      use b2mod_external
      use b2mod_version  ! ajs 07.04.20
#ifdef NAGF95
      use f90_unix       ! IGNORE
#endif
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_clogau
      use eirmod_ctrcei
      use eirmod_ctrig
      use eirmod_cinit
     , , only: eirene_init_cinit, fort_lc
      use eirmod_wneutrals
      use eirmod_infcop
      use eirmod_braeir
      use eirmod_cestim
      use eirmod_coutau
     , , only: eirene_dealloc_coutau
      use eirmod_comusr
     , , only: eirene_dealloc_comusr,
     ,         nmassa, nmassm, nmassi
      use eirmod_comprt
     , , only: iunout
      use eirmod_eirbra
     & , smo_eir => smo
#else
      use b2mod_braeir
      use b2mod_eirbra
     & , smo_eir => smo
#endif
      use b2mod_interp
      use b2mod_subsys
      use b2mod_dimensions
      implicit none
#ifndef NO_CDF
#     include <netcdf.inc>
#endif
      integer nxd,nyd,nsd_local
c
c  version : 27.07.96 20:37
c
      logical ratio
      common /t_ratio/ ratio

      real (kind=R8) ::
     *   pr(-1:nxd,-1:nyd), ro(-1:nxd,-1:nyd,0:nsd_local-1),
     3   energy_balance(-1:nxd,-1:nyd),energy_balance2(-1:nxd,-1:nyd),
     4   pz(-1:nxd,-1:nyd),
     4   energy_balance_old(-1:nxd,-1:nyd),
     3   particle_balance(-1:nxd,-1:nyd,0:nsd_local-1), pb0,
     4   bas_vis_piece(-1:nxd,-1:nyd,0:nsd_local-1),
     *   r_ref, z_ref, s_ref(-1:nyd), profile(-1:nyd),
     7   uedge_c(-1:nxd,-1:nyd), uedge_v(-1:nxd,-1:nyd,1:4),
     8   corners(-1:nxd+1,-1:nyd+1),
     9   zaf(0:nsd_local-1), tif, pxf
      real (kind=R8), allocatable ::
     *   wall_tmp(:,:), trg_tmp(:,:)
      real (kind=R8) ::
     *   fetx(-1:nxd,-1:nyd), fety(-1:nxd,-1:nyd),
     *   fhtx(-1:nxd,-1:nyd), fhty(-1:nxd,-1:nyd),
     *   lnlam(-1:nxd,-1:nyd)
      real (kind=R8) ::
     *   fnax(-1:nxd,-1:nyd,0:nsd_local-1),
     *   fnay(-1:nxd,-1:nyd,0:nsd_local-1),
     *   fmox(-1:nxd,-1:nyd,0:nsd_local-1),
     *   fmoy(-1:nxd,-1:nyd,0:nsd_local-1),
     *   fhmx(-1:nxd,-1:nyd,0:nsd_local-1),
     *   fhmy(-1:nxd,-1:nyd,0:nsd_local-1),
     *   fhpx(-1:nxd,-1:nyd,0:nsd_local-1),
     *   fhpy(-1:nxd,-1:nyd,0:nsd_local-1),
     *   vv(-1:nxd,-1:nyd,0:nsd_local-1),
     *   up(-1:nxd,-1:nyd,0:nsd_local-1),
     *   ww(-1:nxd,-1:nyd,0:nsd_local-1),
     *   safetyfactor(-1:nyd), shear(-1:nyd), rave(-1:nyd),
     1   cvsa_cl_x(-1:nxd,-1:nyd,0:nsd_local-1),
     2   cvsa_cl_y(-1:nxd,-1:nyd,0:nsd_local-1),
     3   div_fna(-1:nxd,-1:nyd,0:nsd_local-1),
     *   div_fmo(-1:nxd,-1:nyd,0:nsd_local-1),
     3   div_fhm(-1:nxd,-1:nyd,0:nsd_local-1),
     *   div_fhp(-1:nxd,-1:nyd,0:nsd_local-1),
     3   div_fhe(-1:nxd,-1:nyd), div_fhi(-1:nxd,-1:nyd),
     3   div_fhet(-1:nxd,-1:nyd), div_fhit(-1:nxd,-1:nyd),
     4   div_fch(-1:nxd,-1:nyd), div_fhj(-1:nxd,-1:nyd),
     5   fllim0fna_x(-1:nxd,-1:nyd,0:nsd_local-1),
     5   fllim0fna_y(-1:nxd,-1:nyd,0:nsd_local-1),
     5   fllim0fhi_x(-1:nxd,-1:nyd,0:nsd_local-1),
     5   fllim0fhi_y(-1:nxd,-1:nyd,0:nsd_local-1)
      logical, save :: first_neoclassical = .true.
      logical, save :: rerun_neoclassical = .true.
      real (kind=R8) ::                                                !cjk
     *   dn_neo(-1:nxd,-1:nyd,0:nsd_local),
     *   chi_neo(-1:nxd,-1:nyd,0:nsd_local),
     *   vn_neo(-1:nxd,-1:nyd,0:nsd_local),
     *   vw_neo(-1:nxd,-1:nyd,0:nsd_local)
      real (kind=R8) ::                                                !cjk
     *   nuab(-1:nxd,-1:nyd,0:nsd_local,0:nsd_local),
     *   nua(-1:nxd,-1:nyd,0:nsd_local),
     *   nustar(-1:nxd,-1:nyd,0:nsd_local),
     *   rhol(-1:nxd,-1:nyd,0:nsd_local),
     *   lnn(-1:nxd,-1:nyd,0:nsd_local),
     *   ltt(-1:nxd,-1:nyd,0:nsd_local),
     *   deltal(-1:nxd,-1:nyd,0:nsd_local),
     *   wb(-1:nxd,-1:nyd,0:nsd_local),
     *   deltaw(-1:nxd,-1:nyd,0:nsd_local),
     *   deltawn(-1:nxd,-1:nyd,0:nsd_local),
     *   deltawt(-1:nxd,-1:nyd,0:nsd_local),
     *   omega(-1:nxd,-1:nyd,0:nsd_local),
     *   deltab(-1:nxd,-1:nyd,0:nsd_local),
     *   deltapar(-1:nxd,-1:nyd,0:nsd_local),
     *   vtha(-1:nxd,-1:nyd,0:nsd_local)
      real (kind=R8), allocatable :: mach(:,:,:)
#ifdef B25_EIRENE
      real (kind=R8), allocatable :: emach(:,:)
#endif
      real (kind=R8) ::
     1   vappr_a, vappr_b, sput_frc, sput_phys, BoRiS,
     2   stick_coef, plate_kb, plate_ka,
     3   q, sput_mixed_alpha, sput_mixed_beta
      real (kind=R8) ::
     *   balance(0:DEF_NATM-1)
      real (kind=R8) ::
     *   mass,vel,vth,press
      integer ispecies(0:nsd_local-1), ismain, ismain0,
     *   neutrals_namelist, sput_frac_flag, use_eirene

#ifndef NO_CDF
      real (kind=R8) ::
     *   movietime(1)
      integer :: iret, frameindex(1), ncidcdf, varid
      integer ntstepcdf, istepcdf
#endif
      integer lucdf
      integer itmp, istart, iend,
     * sput_src, refl_model, sput_phys_model,
     * sput_chem_model, sput_phys_col, drift_style
      character*256 stmp, filename, varname
      character*256 stmp1, stmp2
      integer lu, save_lu(10), ipos, save_ipos(10), lu_pos
      integer is, iss, isp, is2
      character*512 textline, save_textline(10)
      character*8 cval
      character*1 cc
      character*256 sval
      character*256 description
      character*24 hlp_frm
      CHARACTER FILNAM*8, H123*4, REAC*9, CRC*3
c
      real (kind=R8) ::
     *   rshift, zshift, r0, z0, a0, b0, tt0
      common /chord_shift/ rshift, zshift, r0, z0, a0, b0, tt0
c
      logical ex, cflag, found, eirc_done
      logical open_file, streql, isaletter
      integer len_of_digits
      real (kind=R8) ::
     *   qevap, qrad,
     *   Kappa, dKappabydT, Kappa_perp, dKappa_perpbydT, Cp,
     *   mass_density,
     *   layer_Kappa, layer_Cp, layer_emissivity
      external open_file, streql, isaletter,
     *   xertst, xerrab, ipgetr, ipgeti, len_of_digits,
     *   qevap, qrad,
     *   Kappa, dKappabydT, Kappa_perp, dKappa_perpbydT, Cp,
     *   mass_density,
     *   layer_Kappa, layer_Cp, layer_emissivity
      real (kind=R8) ::
     *   t1,t2
      integer it1,it2
      character*256 arg
      integer iflag
#ifdef B25_EIRENE
      character*72 zeile
#endif
c
      real (kind=R8) ::
     * b0r0,dummy,rval,sumjxa,sumjxi,suml(2),sumr(2),time
      real (kind=R8) :: hh, hhn, phi_app
      integer i,j,icell,il,imi,
     1 is1,ival,ix,ix1,ix2,iy,iy1,iy2,izend,izstrt,
     2 ikind,nz,ixl,ixr,iyb,iyt,isel,
     3 i1,ibnd,ixn,iyn,ixf,iyf,idir,isign
      data sput_src/1/, ismain/1/, sput_frc/0.0_R8/,
     5     plate_kb/1.6e-3_R8/, plate_ka/1.3e-5_R8/,
     6     vappr_a/-40181.0_R8/, vappr_b/14.80_R8/, stick_coef/0.4_R8/,
     7     sput_phys_col/3/, sput_chem_model/0/, refl_model/1/,
     8     sput_phys_model/1/, sput_phys/0.0_R8/,
     *     neutrals_namelist/0/, BoRiS/0.0_R8/, sput_frac_flag/0/,
     1     sput_mixed_alpha/1.0_R8/, sput_mixed_beta/1.0_R8/,
     2     use_eirene/0/, ismain0/0/
      save sput_src, sput_frc,
     1     plate_kb, plate_ka, vappr_a, vappr_b, stick_coef, ismain,
     2     refl_model, sput_chem_model,
     3     sput_phys_model, sput_phys_col, sput_phys, BoRiS,
     4     neutrals_namelist,
     5     sput_frac_flag, use_eirene, ismain0, drift_style

c
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
*     ------------------------------------------------------------------

c Function to check the correct parameters!
c
      call subini ('b2plot')
      call ipgeti ('b2mndr_eirene', use_eirene)
#ifdef B25_EIRENE
      CALL EIRENE_DEFAULTS_USR
      filename='input.dat'
      call find_file(filename,ex)
      if (ex.and.use_eirene.ne.0) then
        open(unit=45,file=filename)
        read(45,'(A72)',END=6999) ZEILE
        rewind(45)
        CALL EIRENE_FIND_PARAM(45)
        close(45)
        CALL EIRENE_SET_PARMMOD(1)
        write(*,*) 'Called FIND_PARAM'
        write(*,*) '   NDX = ',NDX
        write(*,*) '   NDY = ',NDY
        write(*,*) ' NSTRA = ',NSTRA
        write(*,*) '   NFL = ',NFL
        write(*,*) '  NATM = ',NATM
        write(*,*) '  NMOL = ',NMOL
        write(*,*) '  NION = ',NION
        write(*,*) 'NLIMPS = ',NLIMPS
      endif
      goto 7000
 6999 continue
      ex = .false.
 7000 continue
      if (.not.ex.or.use_eirene.eq.0) then
        ndx=DEF_NXD+DEF_NCUT*5
        ndy=DEF_NYD
        nfl=DEF_NFL
        natm=DEF_NATM
        nmol=DEF_NMOL
        nion=DEF_NION
        npls=DEF_NPLS
        nlim=DEF_NLIM
        nsts=DEF_NSTS
        nstra=DEF_NSTRA
        nsrfs=DEF_NSRFS
        N1ST=DEF_NYD+1
        N2ND=DEF_NXD+MAX(1,DEF_NCUT/2*3)+
     .       DEF_ISOEXTRA*MAX(0,DEF_NCUT/2-1)
        N3RD=DEF_NPHID
        NREAC=7
        IUNOUT=6
        CALL EIRENE_ALLOC_CLOGAU
        CALL EIRENE_SET_PARMMOD(1)
        CALL EIRENE_INIT_CINIT
      endif
#else
      call xertst(use_eirene.eq.0,
     &    'Compile code with B25_EIRENE option for coupled runs')
#endif
      natmi = natm
      nmoli = nmol
      nioni = nion
      NFLAI = NFL
#ifdef B25_EIRENE
      CALL EIRENE_ALLOC_BRAEIR(NDX,NDY,NFL)
      CALL EIRENE_ALLOC_EIRBRA(NDX,NDY,NFL,NATM,NMOL,NION,NSTRA)
      CALL EIRENE_ALLOC_CTRCEI(1)
      CALL EIRENE_ALLOC_CTRCEI(2)
      CALL EIRENE_SETCON
      if (use_eirene.ne.0) CALL EIRENE_WNEUTRALS_INIT
#else
      call alloc_b2mod_braeir(nxd,nyd,nsd_local)
      call alloc_b2mod_eirbra(nxd,nyd,nsd_local,
     .                        natm,nmol,nion,npls,nstra)
#endif
      call alloc_b2mod_b2plot(nxd,nyd,nsd_local,ns_ext,wklng,
     1 natm,nmol,nion,nlimps)
      allocate(mach(-1:nxd,-1:nyd,
     *     0:max(1,nsd_local,ns_ext,natm,nmol,nion)-1))
      solpstop=' '
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
      call get_environment_variable('SOLPSTOP',
     . status=ierror, length=lenval)
      if (ierror.eq.0) call get_environment_variable('SOLPSTOP',
     . value=solpstop)
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
      texttargsp=''
      do i=1,nparameters
        values(i)=0.0_R8
      end do
      time=tim
#ifndef NO_CDF
      ntstepcdf=0
      istepcdf=0
#endif
      lucdf=0
      eirc_done=.false.
c
      ns=max(nsd_local,1)
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',ns,'i3)'
      write(*,hlp_frm) 'zn: ',(nint(zn(is)),is=0,ns-1)
      do is=0,ns-1
        call species(is,textpl(is),.true.)
        call species(is,textplel(is),.true.)
        call species(is,textwl(is),.true.)
      end do
      textplel(ns)=' el'
      iss=0
      il=0
      do is=0,ns-1
        if(nint(zn(is)) .ne. il) then
          il=nint(zn(is))
          imi=nint(am(is))
          if(il .eq. 0) then
            textsp(iss)=' '//'n'
          else if(il .le. 1) then
            if(imi .eq. 1) textsp(iss)=' '//'H'
            if(imi .eq. 2) textsp(iss)=' '//'D'
            if(imi .eq. 3) textsp(iss)=' '//'T'
          else
            if(imi.ne.nuclear_mass(il)) then
              if(imi.lt.10) then
                write(textsp(iss),'(a,i1,a)') ' ',imi,trim(elements(il))
              else if(imi.lt.100) then
                write(textsp(iss),'(a,i2,a)') ' ',imi,trim(elements(il))
              else
                write(textsp(iss),'(a,i3,a)') ' ',imi,trim(elements(il))
              end if
            else
              textsp(iss)=' '//elements(il)
            end if
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
      textni(0)='All heavies'
      textni(1)='Ions only'
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
      texttrack(0)=' Bulk'

!!!      if(safetyfactor(-1) .eq. 0.0_R8) call init_q

cxpb  This call is moved up from 'f31' so that we used the original cell-faced
cxpb  quantities before their mapping to cell centers
      if (redef_gmtry.eq.0) then
        drift_style = 1
      else
        drift_style = 2
      end if
      call ipgeti ('b2tfnb_drift_style', drift_style)
      call ipgetr ('b2news_BoRiS', BoRiS)
      call b2xpni (nx, ny, ns, na, ni)
      nflai=0
      do is=0,ns-1
        if(is_neutral(is)) nflai=nflai+1
      enddo
      if(nflai.gt.NFL) then
        write(*,*) 'copy_background not called because nflai.gt.NFL',
     1   nflai,NFL
      elseif(nx.gt.NXDD.or.ny.gt.NYDD) then
        write(*,*)
     1   'copy_background not called because NX>NXDD or NY>NYDD',
     1   nx,nxdd,ny,nydd
      else
        call copy_background(nx, ny, ns, BoRiS)
      endif
c
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
      if (plate_option.eq.2) plate_option=3 ! we do not want to recompute the values!
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
       found = .false.
       is = ismain
       do while (is.ge.0 .and. .not.found)
         if (is_neutral(is) .and. nint(zn(is)).eq.nint(zn(ismain))
     &                      .and. nint(am(is)).eq.nint(am(ismain))) then
           ismain0 = is
           found = .true.
         end if
         is = is - 1
       end do
       if (.not.found.and.ismain.ne.1) ismain0 = ismain
       call ipgeti ('b2mwti_ismain0', ismain0)
       call xertst (0.le.ismain0.and.ismain0.lt.ns,
     .   'invalid main neutral species index ismain0')
       call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     .   'invalid main neutral species index ismain0')
      call get_boundary_namelist
      call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)

      call initialise_potentials(ns)
      isp=0
      ispecies(0)=isp
      do is=1,ns-1
        if(nint(zn(is)) .ne. nint(zn(is-1)) .or.
     &     nint(am(is)) .ne. nint(am(is-1))) then
          isp=isp+1
        endif
        ispecies(is)=isp
      enddo
      do is=0,ns-1
        write(*,'(i3,1x,i3,1x,f5.1,1x,f5.1,1x,f5.1,4f11.3)')
     1   is,ispecies(is),zamin(is),zamax(is),zn(is),
     2   potimin(is),potimax(is),potmin(is),potmax(is)
      end do
      call alloc_b2mod_elements(ns, zn, am)
c calculate additional fluxes (fhj, fhm, fhp, fht)
!ajs 07.04.20 only recalculate in older version, now in b2mod_plasma
      if (b2fplasma_version.lt.'03.000.007') then
        call b2xpen (nx, ny, ns,
     1    fna, fna, fne, fch, fhe, fhi, kinrgy, ua, cvsa, rpt,
     2    te, ti, po, BoRiS, fhm, fhp, fhj, fht)
      endif
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
c boundary conditions (species independent)
     1         +b2stbc_she(ix,iy)+b2stbc_shi(ix,iy)
     2         +b2stbc_sch(ix,iy)*po(ix,iy)
c recycling (species independent)
     1         +b2stbr_she(ix,iy)+b2stbr_shi(ix,iy)
     2         +b2stbr_sch(ix,iy)*po(ix,iy)
c other sources (species independent)
     1         +b2stbm_she(ix,iy)+b2stbm_shi(ix,iy)
     2         +b2stbm_sch(ix,iy)*po(ix,iy)
c external sources (species independent)
     1         +ext_she(ix,iy)+ext_shi(ix,iy)
     2         +ext_sch(ix,iy)*po(ix,iy)
            do is=0,ns-1
              energy_balance2(ix,iy)=energy_balance2(ix,iy)
c radiation (electron cooling rate)
     1         -rqahe(ix,iy,is)
c boundary conditions (species dependent)
     3         +b2stbc_smo(ix,iy,is)*ua(ix,iy,is)*(1.0_R8-BoRiS)
     4         -b2stbc_sna(ix,iy,is)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
c recycling (species dependent)
     3         +b2stbr_smo(ix,iy,is)*ua(ix,iy,is)*(1.0_R8-BoRiS)
     4         -b2stbr_sna(ix,iy,is)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
c other sources (species dependent)
     3         +b2stbm_smo(ix,iy,is)*ua(ix,iy,is)*(1.0_R8-BoRiS)
     4         -b2stbm_sna(ix,iy,is)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
c external sources (species dependent)
     3         +ext_smo(ix,iy,is)*ua(ix,iy,is)*(1.0_R8-BoRiS)
     4         -ext_sna(ix,iy,is)*kinrgy(ix,iy,is)*(1.0_R8-BoRiS)
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
     1         b2stbc_sna(ix,iy,is)+b2stbr_sna(ix,iy,is)+
     2         b2stbm_sna(ix,iy,is)+ext_sna(ix,iy,is)
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
      call alloc_b2mod_boundary(nx, ny, ns, nbcd, nnisod)
      if (boundary_namelist.eq.1)
     > call read_b2mod_boundary_namelist(nx, ny, ns, .false.)
      call xertst(nbc.le.nbcd,
     . 'increase size of DEF_NBC in b2mod_dimensions')
      call xertst(nniso.le.nnisod,
     . 'increase size in DEF_NISO in b2mod_dimensions')
      call alloc_b2mod_neutrals
     . (nx, ny, ns, natm, nmol, nion, nstraid, nsts, sput_dst)
      if (neutrals_namelist.eq.1 .or. use_eirene.ne.0)
     > call read_neutrals_namelist(nx, ny, ns, .false.)
      call xertst(nstrat.le.nstraid,
     . 'Increase DEF_NSTRA in b2mod_dimensions')
      if(.not.allocated(iz))
     > call alloc_b2mod_diag(nx,ny,ns,nnreg(0),nstra,natm)
      call alloc_b2mod_user(nx,ns,nlim,nmol,ntns)
cxpb  Bring in wall data and compute oft-needed wall quantities
      write(*,*) 'Calling alloc_b2mod_sputter'
      call alloc_b2mod_sputter(nx,ny,ns)
      call alloc_b2mod_layer(nx,ny,nncut,nlim,nsts)
      call alloc_b2mod_wall_init(nx,ny,nncut)
      call alloc_b2mod_wall(ns,natm,nmol,nion)
      call init_wall(nx,ny)
      call read_b2mod_wall_namelist(99,'    ',found)
      allocate(wall_tmp(1:nwall,0:nsd_local-1),
     *   trg_tmp(1:nwall,1:ndepth_nml))
      call alloc_b2mod_b2plot_wall(nsd_local,natmi,nmoli,nioni,
     *   nwall,ndepth_nml,wklng)
      call alloc_b2mod_b2plot_ank(nwall)
      call alloc_b2plot_wall_loading(nlim,nsgmx)
      call xertst (ndepth_nml.gt.1.or.plate_model.eq.0,
     .  'Inconsistent values of ndepth_nml and plate_model')
      do i = 1, ntargsp
        texttrack(i)=' Deposited '//trim(track_species(i))
        texttargsp(i)=' '//trim(track_species(i))
      end do
      call get_wall_data(nx,ny,ns,BoRiS,found)
cxpb  We now use the data already computed in init_wall
      do ibnd = 1, nbnd
        ix = imapx(ibnd)
        iy = imapy(ibnd)
        Phi_app = get_phi_applied(boundary_namelist,ix,iy)
cxpb  North boundary
        if (topiy(ix,iy).eq.ny+1) then
          idir = 1
          isign = 1
          ixf = ix
          iyf = iy
          ixn = bottomix(ix,iy)
          iyn = bottomiy(ix,iy)
          hh = hy1(ix,iy)
          hhn = hy1(ixn,iyn)
cxpb  East boundary
        else if (rightix(ix,iy).eq.nx+1) then
          idir = 0
          isign = 1
          ixf = ix
          iyf = iy
          ixn = leftix(ix,iy)
          iyn = leftiy(ix,iy)
          hh = hx(ix,iy)
          hhn = hx(ixn,iyn)
cxpb  South Boundary
        elseif (bottomiy(ix,iy).eq.-2) then
          idir = 1
          isign = -1
          ixf = topix(ix,iy)
          iyf = topiy(ix,iy)
          ixn = topix(ix,iy)
          iyn = topiy(ix,iy)
          hh = hy1(ix,iy)
          hhn = hy1(ixn,iyn)
cxpb  Western boundary
        elseif (leftix(ix,iy).eq.-2) then
          idir = 0
          isign = -1
          ixf = rightix(ix,iy)
          iyf = rightiy(ix,iy)
          ixn = rightix(ix,iy)
          iyn = rightiy(ix,iy)
          hh = hx(ix,iy)
          hhn = hx(ixn,iyn)
        endif
        wall_heatflux(ibnd)=max(0.0_R8,
     .    isign*(fht(ixf,iyf,idir)+fhi_ext(ixf,iyf,idir)))/
     .                                 gs(ixf,iyf,idir)
*    ..compute cell-faced quantities
        do is = 0, ns-1
          zaf(is) = (rza(ix,iy,is)*hhn+rza(ixn,iyn,is)*hh)/(hhn+hh)
        enddo
        tif = (ti(ix,iy)*hhn+ti(ixn,iyn)*hh)/(hhn+hh)
        do is = 0, ns-1
          incident_energy(ibnd,is)=
     .     max(0.0_R8,tif/ev+(po(ix,iy)-Phi_app)*zaf(is)) ! xpb: Incident energy (eV)
          wall_flux(ibnd,is)=max(0.0_R8,isign*fna(ixf,iyf,idir,is))
        enddo
        do is = 0, ns_ext-1
          pxf = (pt_ext(ix,iy,is)*hhn+pt_ext(ixn,iyn,is)*hh)/
     &          (hhn+hh)
          wall_heatflux(ibnd)=wall_heatflux(ibnd)+
     &     max(0.0_R8,isign*fa_ext(ixf,iyf,idir,is))/
     &         gs(ixf,iyf,idir)*
     &         ((0.5_R8*mp*am_ext(is)*ua_ext(ix,iy,is)**2+
     &           ta_ext(ix,iy,is))*(1.0_R8-BoRiS) + pxf*ev)
        enddo
      enddo
      ndata = nbnd
      do ibnd = 1, nbnd
        wall_area(ibnd) = plate_area(ibnd)
        wall_temp(ibnd) = target_temp(ibnd,1)
      enddo
      target_stack(1:nbnd,1:ndepth_nml,1) =
     1  target_temp(1:nbnd,1:ndepth_nml)
      target_label(1)='Target temperature (K)'
      call set_targetbarlabel(target_label(1))
c
      call mapx(nx,ny,1,gs(:,:,0),work)
      call mapy(nx,ny,1,gs(:,:,1),work)
      call mapx(nx,ny,1,pbs(:,:,0),work)
      call mapy(nx,ny,1,pbs(:,:,1),work)
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
        call mapx(nx,ny,1,fna_mdf(:,:,0,is),work)
        call mapy(nx,ny,1,fna_mdf(:,:,1,is),work)
        call mapx(nx,ny,1,fnaPSch(:,:,0,is),work)
        call mapy(nx,ny,1,fnaPSch(:,:,1,is),work)
        call mapx(nx,ny,1,fna_fcor(:,:,0,is),work)
        call mapy(nx,ny,1,fna_fcor(:,:,1,is),work)
        call mapx(nx,ny,1,fna_he(:,:,0,is),work)
        call mapy(nx,ny,1,fna_he(:,:,1,is),work)
        if (drift_style.gt.0) then
          call mapx(nx,ny,1,uadia(:,:,0,is),work)
          call mapy(nx,ny,1,uadia(:,:,1,is),work)
          call mapx(nx,ny,1,vadia(:,:,0,is),work)
          call mapy(nx,ny,1,vadia(:,:,1,is),work)
          call mapx(nx,ny,1,wadia(:,:,0,is),work)
          call mapy(nx,ny,1,wadia(:,:,1,is),work)
          call mapx(nx,ny,1,vaecrb(:,:,0,is),work)
          call mapy(nx,ny,1,vaecrb(:,:,1,is),work)
        end if
      enddo
      do is=0,ns_ext-1
        call mapx(nx,ny,1,fa_ext(:,:,0,is),work)
        call mapy(nx,ny,1,fa_ext(:,:,1,is),work)
      enddo
      call mapx(nx,ny,1,fhe(:,:,0),work)
      call mapy(nx,ny,1,fhe(:,:,1),work)
      call mapx(nx,ny,1,fhi(:,:,0),work)
      call mapy(nx,ny,1,fhi(:,:,1),work)
      call mapx(nx,ny,1,fhi_ext(:,:,0),work)
      call mapy(nx,ny,1,fhi_ext(:,:,1),work)
      call mapx(nx,ny,1,fhet(:,:,0),work)
      call mapy(nx,ny,1,fhet(:,:,1),work)
      call mapx(nx,ny,1,fhit(:,:,0),work)
      call mapy(nx,ny,1,fhit(:,:,1),work)
      call mapx(nx,ny,1,fhj(:,:,0),work)
      call mapy(nx,ny,1,fhj(:,:,1),work)
      call mapx(nx,ny,1,fht(:,:,0),work)
      call mapy(nx,ny,1,fht(:,:,1),work)
      call mapx(nx,ny,1,fch(:,:,0),work)
      call mapy(nx,ny,1,fch(:,:,1),work)
      call mapx(nx,ny,1,fchdia(:,:,0),work)      ! cd
      call mapy(nx,ny,1,fchdia(:,:,1),work)      ! dd
      call mapx(nx,ny,1,fch_32(:,:,0),work)
      call mapy(nx,ny,1,fch_32(:,:,1),work)
      call mapx(nx,ny,1,fch_52(:,:,0),work)
      call mapy(nx,ny,1,fch_52(:,:,1),work)
      call mapx(nx,ny,1,fne(:,:,0),work)
      call mapy(nx,ny,1,fne(:,:,1),work)
      call mapx(nx,ny,1,fne_32(:,:,0),work)
      call mapy(nx,ny,1,fne_32(:,:,1),work)
      call mapx(nx,ny,1,fne_52(:,:,0),work)
      call mapy(nx,ny,1,fne_52(:,:,1),work)
      call mapx(nx,ny,1,fni_32(:,:,0),work)
      call mapy(nx,ny,1,fni_32(:,:,1),work)
      call mapx(nx,ny,1,fni_52(:,:,0),work)
      call mapy(nx,ny,1,fni_52(:,:,1),work)
      call mapx(nx,ny,1,fchvispar(:,:,0),work)
      call mapy(nx,ny,1,fchvispar(:,:,1),work)
      call mapx(nx,ny,1,fchvisper(:,:,0),work)
      call mapy(nx,ny,1,fchvisper(:,:,1),work)
      call mapx(nx,ny,1,fchin(:,:,0),work)
      call mapy(nx,ny,1,fchin(:,:,1),work)
      call mapx(nx,ny,1,fch_p(:,:,0),work)           !srv 10.06.09
      call mapy(nx,ny,1,fch_p(:,:,1),work)           !srv 10.06.09
      call mapx(nx,ny,1,fhe_mdf(:,:,0),work)
      call mapx(nx,ny,1,fhe_mdf(:,:,1),work)
      call mapx(nx,ny,1,fhi_mdf(:,:,0),work)
      call mapx(nx,ny,1,fhi_mdf(:,:,1),work)
      call mapx(nx,ny,1,fhePSch(:,:,0),work)
      call mapx(nx,ny,1,fhePSch(:,:,1),work)
      call mapx(nx,ny,1,fhiPSch(:,:,0),work)
      call mapx(nx,ny,1,fhiPSch(:,:,1),work)
      call mapx(nx,ny,1,fchvisq(:,:,0),work)
      call mapx(nx,ny,1,fchvisq(:,:,1),work)
      call mapx(nx,ny,1,fchinert(:,:,0),work)
      call mapx(nx,ny,1,fchinert(:,:,1),work)
      if (drift_style.gt.0) then
        call mapx(nx,ny,1,veecrb(:,:,0),work)
        call mapx(nx,ny,1,veecrb(:,:,1),work)
        call mapx(nx,ny,1,vedia(:,:,0),work)
        call mapx(nx,ny,1,vedia(:,:,1),work)
      end if
      call mapx(nx,ny,1,f_luc_ke,work)
      call mapx(nx,ny,1,f_luc_ki,work)
      call mapx(nx,ny,1,f_luc_et,work)
      call mapx(nx,ny,1,f_luc_sg,work)
      call mapx(nx,ny,1,f_luc_al,work)
      call mapx(nx,ny,1,fllim_ke,work)
      call mapx(nx,ny,1,fllim_ki,work)
      call mapx(nx,ny,1,fllim_al,work)
c
#ifdef B25_EIRENE
      filename=fort_lc//'44'
      call find_file(filename,ex)
#else
      ex = .false.
#endif
      if (.not.ex.or.use_eirene.eq.0) then
        natmi=1
        textan(0)=' ZERO'
        nmoli=1
        textmn(0)=' ZERO'
        nioni=1
        textin(0)=' ZERO'
        nsw=ns
        call eirene_mc_init(nx,ny)
        call alloc_b2mod_diag(nx,ny,ns,nnreg(0),nstraid,natm)
        nlimi=0
        p1=0.0_R8
        p2=0.0_R8
      else
        open(unit=44,file=filename)
        call ngread(ismain,iflag)
        if (natmi.gt.0)
     .   textwl(ns:ns+natmi-1)=textan(0:natmi-1)
        if (nmoli.gt.0)
     .   textwl(ns+natmi:ns+natmi+nmoli-1)=textmn(0:nmoli-1)
        if (nioni.gt.0)
     .   textwl(ns+natmi+nmoli:ns+natmi+nmoli+nioni-1)=textin(0:nioni-1)
        nsw=ns+natmi+nmoli+nioni
      endif

#ifdef B25_EIRENE
      filename=fort_lc//'46'
      call find_file(filename,ex)
      if(ex.and.use_eirene.ne.0) then
        open(unit=46,file=filename)
        call ntread
        allocate(emach(max(1,natmi,nmoli,nioni),ntrii))
        DO I = 1, ntrii
          IF (NATMI.gt.0) THEN
            DO IS = 1, NATMI
              IF (pdena(IS,I).gt.0.0_R8) then
                vxa(IS,I) = vxdena(IS,I)/pdena(IS,I)/
     .                      (nmassa(IS)*mp*1000.0_R8)/100.0_R8
                vya(IS,I) = vydena(IS,I)/pdena(IS,I)/
     .                      (nmassa(IS)*mp*1000.0_R8)/100.0_R8
                vza(IS,I) = vzdena(IS,I)/pdena(IS,I)/
     .                      (nmassa(IS)*mp*1000.0_R8)/100.0_R8
              ELSE
                vxa(IS,I) = 0.0_R8
                vya(IS,I) = 0.0_R8
                vza(IS,I) = 0.0_R8
              ENDIF
              earrowx(I,IS-1) = vxa(IS,I)
              earrowy(I,IS-1) = vya(IS,I)
            END DO
          ELSE
            vxa(1,I) = 0.0_R8
            vya(1,I) = 0.0_R8
            vza(1,I) = 0.0_R8
            earrowx(I,0) = 0.0_R8
            earrowy(I,0) = 0.0_R8
          END IF
          IF (NMOLI.GT.0) THEN
            DO IS = 1, NMOLI
              IF (pdenm(IS,I).gt.0.0_R8) then
                vxm(IS,I) = vxdenm(IS,I)/pdenm(IS,I)/
     .                      (nmassm(IS)*mp*1000.0_R8)/100.0_R8
                vym(IS,I) = vydenm(IS,I)/pdenm(IS,I)/
     .                      (nmassm(IS)*mp*1000.0_R8)/100.0_R8
                vzm(IS,I) = vzdenm(IS,I)/pdenm(IS,I)/
     .                      (nmassm(IS)*mp*1000.0_R8)/100.0_R8
              ELSE
                vxm(IS,I) = 0.0_R8
                vym(IS,I) = 0.0_R8
                vzm(IS,I) = 0.0_R8
              ENDIF
            END DO
          ELSE
            vxm(1,I) = 0.0_R8
            vym(1,I) = 0.0_R8
            vzm(1,I) = 0.0_R8
          END IF
          IF (NIONI.GT.0) THEN
            DO IS = 1, NIONI
              IF (pdeni(IS,I).gt.0.0_R8) then
                vxi(IS,I) = vxdeni(IS,I)/pdeni(IS,I)/
     .                      (nmassi(IS)*mp*1000.0_R8)/100.0_R8
                vyi(IS,I) = vydeni(IS,I)/pdeni(IS,I)/
     .                      (nmassi(IS)*mp*1000.0_R8)/100.0_R8
                vzi(IS,I) = vzdeni(IS,I)/pdeni(IS,I)/
     .                      (nmassi(IS)*mp*1000.0_R8)/100.0_R8
              ELSE
                vxi(IS,I) = 0.0_R8
                vyi(IS,I) = 0.0_R8
                vzi(IS,I) = 0.0_R8
              ENDIF
            END DO
          ELSE
            vxi(1,I) = 0.0_R8
            vyi(1,I) = 0.0_R8
            vzi(1,I) = 0.0_R8
          END IF
        END DO
      end if
#else
      ex = .false.
#endif
c
      nh=0.0_R8
      do is=0,ns-1
        stmp1=''
        i=1
        j=0
        do while (i.le.len(textpl(is)))
          cc=textpl(is)(i:i)
          if (cc.eq.'^'.or.cc.eq.'{'.or.cc.eq.'}'.or.
     &       (cc.eq.' '.and.i.gt.1) ) then
            i=i+1
          else if (cc.ge.'A'.and.cc.le.'Z') then
            i=i+1
            j=j+1
            stmp1(j:j)=char(ichar(cc)-ichar('A')+ichar('a'))
          else
            i=i+1
            j=j+1
            stmp1(j:j)=cc
          endif
        enddo
        if(is_neutral(is) .and. nint(zn(is)).eq.1) then
          nh(-1:nx,-1:ny)=nh(-1:nx,-1:ny)+na(-1:nx,-1:ny,is)
          do iss=1,natmi
            stmp2=trim(textan(iss-1))//'0'
            call lower_case(stmp2)
            if(streql(stmp1,stmp2)) then
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
        XX(IX,2) = real(ix,R8)
        XX(IX,1) = real(IX,R8)-0.5_R8
      end do
      do iy=-1,ny
        YY(IY,2) = real(IY,R8)
        YY(IY,1) = real(IY,R8)-0.5_R8
      end do
      do ibnd=1,nbnd
        targetx(ibnd,2) = real(ibnd,R8)
        targetx(ibnd,1) = real(ibnd,R8)-0.5_R8
      end do
      do i=1,ndepth_nml
        targety(i,2) = real(i,R8)
        targety(i,1) = real(i,R8)-0.5_R8
      end do

      call b2xppr(nx,ny,ns,ne,na,te,ti,pr)
      call b2xppz(nx,ny,ns,ne,na,te,ti,pz)
      call b2xzef(nx,ny,ns,rz2,na,ne,zeff)
      psi_sep=0.5_R8*(fpsi(jxa,jsep,2)+fpsi(jxa,jsep,3))
      psi=0.25_R8*sum(fpsi,3)-psi_sep
      fbz=0.25_R8*sum(ffbz,3)
      targetr=0.5_R8*sum(targetcr,3)
      targetz=0.5_R8*sum(targetcz,3)
      icell=0
      r0=0.0_R8
      z0=0.0_R8
      do ix = -1, nx
        if (on_closed_surface(ix,-1)) then
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
      values(1)=minval(xx)                               !xmin
      values(2)=maxval(xx)+0.5_R8                        !xmax
      values(3)=minval(yy)                               !ymin
      values(4)=maxval(yy)+0.5_R8                        !ymax
#ifdef B25_EIRENE
      if (ex.and.use_eirene.ne.0) then
        values(5)=min(minval(crx),minval(xtrian))        !rmin
        values(6)=max(maxval(crx),maxval(xtrian))        !rmax
        values(7)=min(minval(cry),minval(ytrian))        !zmin
        values(8)=max(maxval(cry),maxval(ytrian))        !zmax
      else
#endif
        values(5)=minval(crx)                            !rmin
        values(6)=maxval(crx)                            !rmax
        values(7)=minval(cry)                            !zmin
        values(8)=maxval(cry)                            !zmax
#ifdef B25_EIRENE
      end if
#endif
      values(11)=init_value                              !omin
      values(12)=init_value                              !omax
      values(13)=init_value                              !xref
      values(16)=5.0e-4_R8                               !traj: plottime
      values(17)=0.0_R8                                  !rshift
      values(18)=0.0_R8                                  !zshift
      values(19)=0.0_R8                                  !r0
      values(20)=0.0_R8                                  !z0
      values(21)=0.0_R8                                  !a0
      values(22)=0.0_R8                                  !b0
      values(23)=0.0_R8                                  !tt0
      values(24)=minval(targetx(1:nbnd,:))               !xmin target
      values(25)=maxval(targetx(1:nbnd,:))+0.5_R8        !xmax target
      values(26)=minval(targety(1:ndepth_nml,:))         !ymin target
      values(27)=maxval(targety(1:ndepth_nml,:))+0.5_R8  !ymax target
      values(28)=minval(targetcr(1:nbnd,1:ndepth_nml,:)) !rmin target
      values(29)=maxval(targetcr(1:nbnd,1:ndepth_nml,:)) !rmax target
      values(30)=minval(targetcz(1:nbnd,1:ndepth_nml,:)) !zmin target
      values(31)=maxval(targetcz(1:nbnd,1:ndepth_nml,:)) !zmax target
      values(32)=minval(target_temp(1:nbnd,1:ndepth_nml))!fmin target
      values(33)=maxval(target_temp(1:nbnd,1:ndepth_nml))!fmax target
      values(34)=init_value                              !omin target
      values(35)=init_value                              !omax target
      values(36)=init_value                              !smin wall
      values(37)=init_value                              !smax wall
      values(38)=init_value                              !Eirene min
      values(39)=init_value                              !Eirene max
      expfile=' '
      wallfile=' '
c
      call b2page(1,1)
      color_set=.true.
      call contourplot(.true.)     ! draw contour plot by default
      call nofill(.false.)         ! fill contours by default
      call ipolate(.false.)        ! do not use splines unless wanted
      call clvl(.true.)            ! contour level lines plotted
      call labelon(.true.)         ! contour labels drawn
      call gridon(.false.)         ! overlayed mesh off by default
      call trigon(.false.)         ! overlayed triangular grid off by default
      call vesselon(.false.)       ! same for vessel
      call separatrixon(.false.)   ! same for separatrix
      call adjust_axes(.false.)    ! nice axes by default
      call adjust_bold(.false.)    ! do not bold axes by default
      call set_full_psi_plotting(.false.) ! do not plot against full psi
      call set_norm_psi_plotting(.false.) ! do not plot against normalised psi
      call set_root_psi_plotting(.false.) ! do not plot against root(psi)
      call set_midplane_mapping(.false.)  ! do not map to midplane by default
      call wallon(.false.,0,values(14),values(15),1)
      call targeton(.false.,values(32),values(33),1)
      call set_partition_colors(0.0_R8)
      plotting_sep=.false.
      plotting_wall=.false.
      b2plot_debug=0
      sepwidth=1.0_R8            ! standard separatrix line width
      vwidth=1.0_R8              ! standard vessel line width
      owidth=1.0_R8              ! standard outline line width
      ewidth=1.0_R8              ! standard experimental data line width
      arroww=1.0_R8              ! standard arrow line width
      lwidth=1.0_R8              ! linewidth
      wwidth=10.0_R8             ! wall linewidth
      zwidth=1.0_R8              ! axis linewidth
      lfsz=0                     ! line label font size
      typeofplot=3               ! physical and computational
      xparallel=.false.          ! default to f.y using poloidal distance
      expline=.true.             ! default to plotting experimental data as lines
      expsymbol=5                ! default symbol is crosses
      call set_fill_option(0)
      linlog=1                   ! linear plot
      wlinlog=1                  ! linear plot
      tlinlog=0                  ! log plot for trajectories
      trglinlog=1                ! linear plot for targets
      call set_square            ! square plot
      writing_to_file=.false.
      call initialise_plotting_0
      call linewidth(lwidth)
      call axislinewidth(zwidth)
      call set_line_font_size(lfsz)
      moviename=' '
      movie_skip=1
      average=1
      use_compile=.true.
      acol=1
      bcol=1
      dcol=1
      ecol=1
      fcol=1
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
      neo_ic=3    ! banana plateau + Pfirsch-Schlueter
      call ipgeti ('neoclassical_ic', neo_ic)
      call xertst (0.le.neo_ic .and. neo_ic.le.4,
     . 'faulty parameter neoclassical_ic')
      set_resignore_regions_to_zero=.false.
      ue(:,:)=-fch_p(:,:,0)/pbs(:,:,0)/qe       !srv 10.06.09
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            ro(ix,iy,is)=am(is)*mp*na(ix,iy,is)
            vv(ix,iy,is)=fna(ix,iy,1,is)/na(ix,iy,is)/gs(ix,iy,1)
            up(ix,iy,is)=
     .       fna(ix,iy,0,is)/na(ix,iy,is)/gs(ix,iy,0)/qz(ix,iy,1)
            ue(ix,iy)=ue(ix,iy)+rza(ix,iy,is)*ua(ix,iy,is)*na(ix,iy,is)
            fnax(ix,iy,is)=fna(ix,iy,0,is)
            fmox(ix,iy,is)=fmo(ix,iy,0,is)
            fhmx(ix,iy,is)=fhm(ix,iy,0,is)
            fhpx(ix,iy,is)=fhp(ix,iy,0,is)
            arrowx(ix,iy,is)=up(ix,iy,is)
            fnay(ix,iy,is)=fna(ix,iy,1,is)
            fmoy(ix,iy,is)=fmo(ix,iy,1,is)
            fhmy(ix,iy,is)=fhm(ix,iy,1,is)
            fhpy(ix,iy,is)=fhp(ix,iy,1,is)
            arrowy(ix,iy,is)=vv(ix,iy,is)
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
      ue(:,:)=ue(:,:)+ue_ext(:,:)*ne_ext(:,:)
      ue(:,:)=ue(:,:)/ne(:,:)
      do iy=-1,ny
        rave(iy)=0.0_R8
        t1=0.0_R8
        t2=0.0_R8
        icell=0
        q=0.0_R8
        do ix=-1,nx
          if(on_closed_surface(ix,iy)) then
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
        if (safetyfactor(iy).eq.0.0_R8 .or. iy.gt.jsep) then
          shear(iy) = 0.0_R8
        else if (iy.eq.-1) then
          shear(iy)=(safetyfactor(iy+1)-safetyfactor(iy))/
     &     (rave(iy+1)-rave(iy))
        elseif (iy.lt.jsep) then
          shear(iy)=(safetyfactor(iy+1)-safetyfactor(iy-1))/
     &     (rave(iy+1)-rave(iy-1))
        elseif (iy.eq.jsep) then
          shear(iy)=(safetyfactor(iy)-safetyfactor(iy-1))/
     &     (rave(iy)-rave(iy-1))
        endif
      enddo
c
      lu=5
      lu_pos=0
      call readline(textline,lu,ipos)
c
      irstack=0
      iistack=0
      isstack=0
      iwrk=0
      iewrk=0
      iwwrk=0
      itwrk=1

   10 continue
      if(b2plot_debug.ge.5)
     1 write(*,'(a,7i4)') 'STACK i,r,s,m,e,w,t = ',
     . iistack,irstack,isstack,iwrk,iewrk,iwwrk,itwrk
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
        if(b2plot_debug.ge.6) then
          if(streql(sval,' ')) then
            write(*,*) 'Read empty string'
          else
            write(*,*) 'Read string: ',trim(sval)
          endif
        endif
        goto 10
      end if
      call lower_case(cval)                                !srv 14.10.13

      cflag=.false.
      if (.not.cflag) call select_matrix
      if (.not.cflag) call select_stack_manipulation
      if (.not.cflag) call select_parameters
      if (.not.cflag) call select_command
      if (.not.cflag) call select_eirene_stack
      if (.not.cflag) call select_eirene_command
      if (.not.cflag) call select_wall_stack
      if (.not.cflag) call select_wall_command
      if (.not.cflag) call select_target_stack
      if (.not.cflag) call select_constants
      if (.not.cflag) call select_flvalue

      if (.not.cflag) then
        if (cval.eq.'eof') then    ! end of file
          if (lu.eq.5) then
            if (writing_to_file) close(nwtf)
            call terminate_plotting
            write(*,'(a,7i4)') 'STACK i,r,s,m,e,w,t = ',
     .       iistack,irstack,isstack,iwrk,iewrk,iwwrk,itwrk
            deallocate(mach)
            call dealloc_b2mod_diag
            call dealloc_b2mod_wall
            call dealloc_b2mod_layer
            call dealloc_b2mod_b2plot
            call dealloc_sputter_data
            call dealloc_b2mod_sputter
            call dealloc_b2mod_neutrals
            call dealloc_b2mod_elements
            call dealloc_b2mod_boundary
            call dealloc_b2mod_wall_init
            call dealloc_b2mod_b2plot_xy
            call dealloc_b2mod_b2plot_ank
            call dealloc_b2mod_b2plot_wall
            call dealloc_b2mod_neoclassical
            call dealloc_b2plot_wall_loading
#ifdef B25_EIRENE
            if (allocated(emach)) then
              deallocate(emach)
              call dealloc_b2mod_b2plot_eirene
            end if
            CALL EIRENE_DEALLOC_CTRIG
            CALL EIRENE_DEALLOC_CESTIM
            CALL EIRENE_DEALLOC_COMUSR
            CALL EIRENE_DEALLOC_COUTAU
            CALL EIRENE_DEALLOC_CTRCEI
            CALL EIRENE_DEALLOC_BRAEIR
            CALL EIRENE_DEALLOC_EIRBRA
            CALL EIRENE_DEALLOC_INFCOP
            CALL EIRENE_EXTRAB25_CLEANUP
            CALL EIRENE_WNEUTRALS_DEALLOC
#else
            call dealloc_b2mod_braeir
            call dealloc_b2mod_eirbra
            call dealloc_b2mod_eirdiag
#endif
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
          if (writing_to_file) close(nwtf)
          call terminate_plotting
          write(*,'(a,7i4)') 'STACK i,r,s,m,e,w,t = ',
     .     iistack,irstack,isstack,iwrk,iewrk,iwwrk,itwrk
          deallocate(mach)
          call dealloc_b2mod_diag
          call dealloc_b2mod_wall
          call dealloc_b2mod_layer
          call dealloc_b2mod_b2plot
          call dealloc_sputter_data
          call dealloc_b2mod_sputter
          call dealloc_b2mod_elements
          call dealloc_b2mod_boundary
          call dealloc_b2mod_neutrals
          call dealloc_b2mod_wall_init
          call dealloc_b2mod_b2plot_xy
          call dealloc_b2mod_b2plot_ank
          call dealloc_b2mod_b2plot_wall
          call dealloc_b2mod_neoclassical
          call dealloc_b2plot_wall_loading
#ifdef B25_EIRENE
          if (allocated(emach)) then
            deallocate(emach)
            call dealloc_b2mod_b2plot_eirene
          end if
          CALL EIRENE_DEALLOC_CTRIG
          CALL EIRENE_DEALLOC_CESTIM
          CALL EIRENE_DEALLOC_COMUSR
          CALL EIRENE_DEALLOC_COUTAU
          CALL EIRENE_DEALLOC_CTRCEI
          CALL EIRENE_DEALLOC_BRAEIR
          CALL EIRENE_DEALLOC_EIRBRA
          CALL EIRENE_DEALLOC_INFCOP
          CALL EIRENE_EXTRAB25_CLEANUP
          CALL EIRENE_WNEUTRALS_DEALLOC
#else
          call dealloc_b2mod_braeir
          call dealloc_b2mod_eirbra
          call dealloc_b2mod_eirdiag
#endif
          call subend ()
          return
        else
          write(*,*) 'Unknown ',trim(cval)
        endif
      endif
      goto 10

      contains

      integer function ixx(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      ixx = mod(ix+1,nx+2)-1
      return
      end function ixx

      integer function iyy(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      iyy = iy+(ix+1)/(nx+2)
      return
      end function iyy

      subroutine select_matrix
      integer ix, iy, is, iso, ibc

      cflag=.true.
      CVALUE: select case (cval)
        case('pit','pitch')
          description='Pitch'
          nz=1
          do ix = -1, nx
            do iy = -1, ny
              mach(ix,iy,0)=bb(ix,iy,0)/bb(ix,iy,3)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case('bx','bpol')
          description='b_{x} (T) (poloidal)'
          nz=1
          call copy(bb(:,:,0),1.0_R8,nz,textbl)
        case('by','brad')
          description='b_{y} (T)'
          nz=1
          call copy(bb(:,:,1),1.0_R8,nz,textbl)
        case('bz','btor')
          description='b_{z} (T) (toroidal)'
          nz=1
          call copy(bb(:,:,2),1.0_R8,nz,textbl)
        case('bb')
          description='bb (T)'
          nz=1
          call copy(bb(:,:,3),1.0_R8,nz,textbl)
        case('b')
          description='B (T) (components)'
          nz=4
          call copy(bb,1.0_R8,nz,textcomponents)
        case('ffbz')
          description='Toroidal magnetic flux function (Wb/m)'
          nz=1
          call copy(fbz,1.0_R8,nz,textbl)
        case('psi')
          description='Poloidal magnetic flux function (Wb)'
          nz=1
          call copy(psi,1.0_R8,nz,textbl)
        case('na')
          description='Ion density (m^{-3})'
          nz=ns
          call copy(na,1.0_R8,nz,textpl)
        case('nh')
          description='Hydrogenic atom density (m^{-3})'
          nz=1
          call copy(nh,1.0_R8,nz,textbl)
        case('ni')
          description='Total density (m^{-3})'
          nz=2
          call copy(ni,1.0_R8,nz,textni)
        case('ua')
          description='Parallel velocity (u_{a}) (m.s^{-1})'
          nz=ns
          call copy(ua,1.0_R8,nz,textpl)
        case('ue')
          description='Parallel electron velocity (u_{e}) (m.s^{-1})'
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
                ww(ix,iy,is)= ua(ix,iy,is)*bb(ix,iy,2)/bb(ix,iy,3) -
     .             (vaecrb(ix,iy,0,is)+wadia(ix,iy,0,is))*
     .             bb(ix,iy,0)/bb(ix,iy,2)
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
          description='Static pressure (Pa)'
          nz=1
          call copy(pr,1.0_R8,nz,textbl)
        case('pz')
          description='Plasma pressure (Pa)'
          nz=1
          call copy(pz,1.0_R8,nz,textbl)
        case('ne')
          description='Electron density (m^{-3})'
          nz=1
          call copy(ne,1.0_R8,nz,textbl)
        case('ne2')
          description='ne2 (m^{-3})'
          nz=1
          mach(-1:nx,-1:ny,0)=ne_ext
          do is = 0, ns-1
             mach(-1:nx,-1:ny,0) = mach(-1:nx,-1:ny,0) +
     &            rz2(-1:nx,-1:ny,is)*na(-1:nx,-1:ny,is)
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case('po')
          description='Electric potential (V)'
          nz=1
          call copy(po,1.0_R8,nz,textbl)
        case('dfhe')
          description='\mtD\dfhe (W)'
          nz=1
          call copy(div_fhe,1.0_R8,nz,textbl)
        case('dfhi')
          description='\mtD\dfhi (W)'
          nz=1
          call copy(div_fhi,1.0_R8,nz,textbl)
        case('dhet','dfhet')
          description='\mtD\dfhet (W)'
          nz=1
          call copy(div_fhet,1.0_R8,nz,textbl)
        case('dhit','dfhit')
          description='\mtD\dfhit (W)'
          nz=1
          call copy(div_fhit,1.0_R8,nz,textbl)
        case('dfhj')
          description='\mtD\dfhj (W)'
          nz=1
          call copy(div_fhj,1.0_R8,nz,textbl)
        case('dfch')
          description='\mtD\dfch (A)'
          nz=1
          call copy(div_fch,1.0_R8,nz,textbl)
        case('dfna')
          description='\mtD\dfna (particles/s)'
          nz=ns
          call copy(div_fna,1.0_R8,nz,textpl)
        case('dfmo')
          description='\mtD\dfmo (N.s)'
          nz=ns
          call copy(div_fmo,1.0_R8,nz,textpl)
        case('dfhm')
          description='\mtD\dfhm (W)'
          nz=ns
          call copy(div_fhm,1.0_R8,nz,textpl)
        case('dfhp')
          description='\mtD\dfhp (W)'
          nz=ns
          call copy(div_fhp,1.0_R8,nz,textpl)
        case('zeff')
          description='Z_{eff}'
          nz=1
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
        case('fnex')
          description='Poloidal electron flux (s^{-1})'
          nz=1
          call copy(fne(:,:,0),1.0_R8,1,textbl)
        case('fney')
          description='Radial electron flux (s^{-1})'
          nz=1
          call copy(fne(:,:,1),1.0_R8,1,textbl)
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
          fhtx(:,:)=fht(:,:,0)+fhi_ext(:,:,0)
          do is=0,ns_ext-1
            fhtx(:,:)=fhtx(:,:)+fa_ext(:,:,0,is)*
     .       (pt_ext(:,:,is)*ev+
     .       (0.5_R8*am_ext(is)*mp*ua_ext(:,:,is)**2+ta_ext(:,:,is))*
     .       (1.0_R8-BoRiS))
          enddo
          call copy(fhtx,1.0_R8,nz,textbl)
        case('fhty')
          description='Radial total energy flux (W)'
          nz=1
          fhty(:,:)=fht(:,:,1)+fhi_ext(:,:,1)
          do is=0,ns_ext-1
            fhty(:,:)=fhty(:,:)+fa_ext(:,:,1,is)*
     .       (pt_ext(:,:,is)*ev+
     .       (0.5_R8*am_ext(is)*mp*ua_ext(:,:,is)**2+ta_ext(:,:,is))*
     .       (1.0_R8-BoRiS))
          enddo
          call copy(fhty,1.0_R8,nz,textbl)
        case('fle')
          description='Electron heat flow flux limit factor'
          nz=1
          call copy(fllime,1.0_R8,nz,textbl)
        case('fli')
          description='Ion heat flow flux limit factor'
          nz=1
          call copy(fllimi,1.0_R8,nz,textbl)
        case('nflfnax')
          description='Poloidal neutral particle flux limit factor'
          nz=ns
          call copy(fllim0fna_x,1.0_R8,nz,textpl)
        case('nflfnay')
          description='Radial neutral particle flux limit factor'
          nz=ns
          call copy(fllim0fna_y,1.0_R8,nz,textpl)
        case('nflfhix')
          description='Poloidal neutral ion energy flux limit factor'
          nz=ns
          call copy(fllim0fhi_x,1.0_R8,nz,textpl)
        case('nflfhiy')
          description='Radial neutral ion energy flux limit factor'
          nz=ns
          call copy(fllim0fhi_y,1.0_R8,nz,textpl)
        case('flvisc')
          description='Viscosity flux limit factor'
          nz=ns
          call copy(fllimvisc,1.0_R8,nz,textpl)
        case('fluce')
          description=
     .     'Electron heat diffusivity Luciani flux limit factor'
          nz=1
          call copy(f_luc_ke,1.0_R8,nz,textbl)
        case('fluci')
          description='Ion heat diffusivity Luciani flux limit factor'
          nz=1
          call copy(f_luc_ki,1.0_R8,nz,textbl)
        case('flucv')
          description='Viscosity Luciani flux limit factor'
          nz=1
          call copy(f_luc_et,1.0_R8,nz,textbl)
        case('flucs')
          description='Electric conductivity Luciani flux limit factor'
          nz=1
          call copy(f_luc_sg,1.0_R8,nz,textbl)
        case('fluca')
          description='Thermo-electric Luciani flux limit factor'
          nz=1
          call copy(f_luc_al,1.0_R8,nz,textbl)
        case('flke')
          description=
     .     'Electron heat diffusivity classical flux limit factor'
          nz=1
          call copy(fllim_ke,1.0_R8,nz,textbl)
        case('flki')
          description=
     .     'Ion heat diffusivity classical flux limit factor'
          nz=1
          call copy(fllim_ki,1.0_R8,nz,textbl)
        case('flal')
          description='Thermo-electric classical flux limit factor'
          nz=1
          call copy(fllim_al,1.0_R8,nz,textbl)
        case('flalc')
          description='Thermo-electric coefficient limit factor'
          nz=1
          call copy(fllim_al_c,1.0_R8,nz,textbl)
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
          call copy(energy_balance,1.0_R8,nz,textbl)
       case('ebls')
          description='Energy balance (sources) (W)'
          nz=1
          call copy(energy_balance2,1.0_R8,nz,textbl)
       case('ebal')
          description='Energy balance (W)'
          nz=1
          call copy(energy_balance+energy_balance2,1.0_R8,nz,textbl)
       case('pbal')
          description='Particle balance (particles)'
          nz=ns
          call copy(particle_balance,1.0_R8,nz,textpl)
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
          mach(:,:,0)=gs(:,:,0)*qz(:,:,1)
          call copy(mach,1.0_R8,nz,textbl)
        case('sy')
          description='sy (m^{2})'
          nz=1
          call copy(gs(:,:,1),1.0_R8,nz,textbl)
        case('sz')
          description='sz (m^{2})'
          nz=1
          call copy(gs(:,:,2),1.0_R8,nz,textbl)
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
          call copy(qc,1.0_R8,nz,textbl)
        case('mach')
          description='Mach number'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              vth=0.0_R8
              vel=0.0_R8
              do is=0,ns-1
                if (is_neutral(is)) cycle
                vth=vth+na(ix,iy,is)*
     .           (rza(ix,iy,is)*te(ix,iy)+ti(ix,iy))/
     .           (am(is)*(ni(ix,iy,1)+ni_ext(ix,iy,1))*mp)
                vel=vel+ua(ix,iy,is)*rza(ix,iy,is)*na(ix,iy,is)/
     .            ne(ix,iy)
              end do
              do is=0,ns_ext-1
                if (is_ext_neutral(is)) cycle
                vth=vth+na_ext(ix,iy,is)*
     .           (za_ext(ix,iy,is)*te(ix,iy)+ta_ext(ix,iy,is))/
     .           (am_ext(is)*(ni(ix,iy,1)+ni_ext(ix,iy,1))*mp)
                vel=vel+
     .            ua_ext(ix,iy,is)*za_ext(ix,iy,is)*na_ext(ix,iy,is)/
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
              do is=0,ns_ext-1
                if (.not.is_ext_neutral(is)) then
                  mass=mass+am_ext(is)*mp*na_ext(ix,iy,is)
                  vel=vel+
     .             za_ext(ix,iy,is)*na_ext(ix,iy,is)*ua_ext(ix,iy,is)
                  press=press+na_ext(ix,iy,is)*ta_ext(ix,iy,is)
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
          call b2xehx(nx, ny, te, po, ne, mach(-1,-1,0))
          call copy(mach,1.0_R8,nz,textbl)
        case('ehy')
          description='ehy'
          nz=1
          call b2xehy(nx, ny, te, po, ne, mach(-1,-1,0))
          call copy(mach,1.0_R8,nz,textbl)
        case('dab2')
          description='Atomic density (m^{-3})'
          if (natmi.gt.0) then
            nz=natmi
            mach(:,:,0:nz-1)=dab2(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textan)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case('dmb2')
          description='Molecular density (m^{-3})'
          if (nmoli.gt.0) then
            nz=nmoli
            mach(:,:,0:nz-1)=dmb2(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textmn)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case('ahal')
          description=
     1     'Atomic H-\Ga\R emissivity (photons.m^{-3}.s^{-1})'
          nz=1
          mach(:,:,0)=emiss(0:nx+1,0:ny+1,1,1)
          call copy(mach,1.0_R8,nz,textbl)
        case ('mhal')
          description=
     1     'Molecular H-\Ga\R emissivity (photons.m^{-3}.s^{-1})'
          nz=1
          mach(:,:,0)=emissmol(0:nx+1,0:ny+1,1,1)
          call copy(mach,1.0_R8,nz,textbl)
        case ('srcml')
          description='Molecule dissociation rate (#.s^{-1})'
          if (nmoli.gt.0) then
            nz=nmoli
            mach(:,:,0:nz-1)=srcml(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,-1.0_R8/qe,nz,textmn)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case ('edissml')
          description='Power to atoms due to molecule dissociation (W)'
          if (nmoli.gt.0) then
            nz=nmoli
            mach(:,:,0:nz-1)=edissml(0:nx+1,0:ny+1,1:nz,0)
            call copy(mach,-1.0_R8,nz,textmn)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          endif
        case ('neutrad')
          description='Power loss due to atom radiation (W)'
          if (natmi.gt.0) then
            nz=natmi
            mach(:,:,0:nz-1)=eneutrad(0:nx+1,0:ny+1,1:nz,0)
            call copy(mach,-1.0_R8,nz,textan)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          endif
        case ('molrad')
          description='Power loss due to molecule radiation (W)'
          if (nmoli.gt.0) then
            nz=nmoli
            mach(:,:,0:nz-1)=emolrad(0:nx+1,0:ny+1,1:nz,0)
            call copy(mach,-1.0_R8,nz,textmn)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          endif
        case ('ionrad')
          description='Power loss due to molecular ion radiation (W)'
          if (nioni.gt.0) then
            nz=nioni
            mach(:,:,0:nz-1)=eionrad(0:nx+1,0:ny+1,1:nz,0)
            call copy(mach,-1.0_R8,nz,textin)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          endif
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
     1     'Ionisation parallel momentum losses (N}) [rsamo]'
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
     1     'Recombination parallel momentum losses (N) [rramo]'
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
     1     'Ion parallel momentum neutral losses (N) [rcxmo]'
          nz=ns
          call copy(rcxmo,1.0_R8,nz,textpl)
        case ('fetx')
          description='Total poloidal internal energy flux (W)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              fetx(ix,iy)=fht(ix,iy,0)-fhj(ix,iy,0)+fhi_ext(ix,iy,0)
              do is=0,ns_ext-1
                fetx(ix,iy)=fetx(ix,iy)+
     .           pt_ext(ix,iy,is)*ev*fa_ext(ix,iy,0,is)+
     .          (0.5_R8*am_ext(is)*mp*ua_ext(ix,iy,is)**2+
     .           ta_ext(ix,iy,is))*fa_ext(ix,iy,0,is)*(1.0_R8-BoRiS)
              end do
            end do
          end do
          call copy(fetx,1.0_R8,nz,textbl)
        case ('fety')
          description='Total radial internal energy flux (W)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              fety(ix,iy)=fht(ix,iy,1)-fhj(ix,iy,1)+fhi_ext(ix,iy,1)
              do is=0,ns_ext-1
                fety(ix,iy)=fety(ix,iy)+
     .           pt_ext(ix,iy,is)*ev*fa_ext(ix,iy,1,is)+
     .          (0.5_R8*am_ext(is)*mp*ua_ext(ix,iy,is)**2+
     .           ta_ext(ix,iy,is))*fa_ext(ix,iy,1,is)*(1.0_R8-BoRiS)
              end do
            end do
          end do
          call copy(fety,1.0_R8,nz,textbl)
        case ('tab2')
          description='Neutral atom temperature (eV)'
          if (natmi.gt.0) then
            nz=natmi
            mach(:,:,0:nz-1)=tab2(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8/ev,nz,textan)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case ('tmb2')
          description='Neutral molecule temperature (eV)'
          if (nmoli.gt.0) then
            nz=nmoli
            mach(:,:,0:nz-1)=tmb2(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8/ev,nz,textmn)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case ('rfla')
          description='Radial atomic flux (s^{-1}.m^{-2})'
          if (natmi.gt.0) then
            nz=natmi
            mach(:,:,0:nz-1)=rfluxa(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textan)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case ('rflm')
          description='Radial molecular flux (s^{-1}.m^{-2})'
          if (nmoli.gt.0) then
            nz=nmoli
            mach(:,:,0:nz-1)=rfluxm(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textmn)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case ('refa')
          description='Radial atomic energy flux (W.m^{-2})'
          if (natmi.gt.0) then
            nz=natmi
            mach(:,:,0:nz-1)=refluxa(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textan)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case ('refm')
          description='Radial molecular energy flux (W.m^{-2})'
          if (nmoli.gt.0) then
            nz=nmoli
            mach(:,:,0:nz-1)=refluxm(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textmn)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case ('pfla')
          description='Poloidal atomic flux (s^{-1}.m^{-2})'
          if (natmi.gt.0) then
            nz=natmi
            mach(:,:,0:nz-1)=pfluxa(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textan)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case ('pflm')
          description='Poloidal molecular flux (s^{-1}.m^{-2})'
          if (nmoli.gt.0) then
            nz=nmoli
            mach(:,:,0:nz-1)=pfluxm(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textmn)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case ('pefa')
          description='Poloidal atomic energy flux (W.m^{-2})'
          if (natmi.gt.0) then
            nz=natmi
            mach(:,:,0:nz-1)=pefluxa(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textan)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
        case ('pefm')
          description=
     1             'Poloidal molecular energy flux (W.m^{-2})'
          if (nmoli.gt.0) then
            nz=nmoli
            mach(:,:,0:nz-1)=pefluxm(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textmn)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          end if
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
          description='Mass density (kg.m^{-3})'
          nz=ns
          call copy(ro,1.0_R8,nz,textpl)
        case ('b2espcr')
          description='Nuclear sequence index'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=real(b2espcr(is),R8)
              end do
            end do
          end do
          call copy(mach,1.0_R8,nz,textpl)
        case ('b2eatcr')
          description='Eirene atomic species index'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                mach(ix,iy,is)=real(b2eatcr(is),R8)
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
          description='Ion charge (e)'
          nz=ns
          call copy(rza,1.0_R8,nz,textpl)
        case ('za2')
          description='Squared ion charge (e^{2})'
          nz=ns
          call copy(rz2,1.0_R8,nz,textpl)
        case ('zn')
          description='Nuclear charge (e)'
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
          description='Dynamic plasma pressure (Pa)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mass=0.0_R8
              vel=0.0_R8
              do is=0,ns-1
                if (is_neutral(is)) cycle
                mass=mass+mp*am(is)*na(ix,iy,is)
                vel=vel+rza(ix,iy,is)*na(ix,iy,is)*ua(ix,iy,is)
              end do
              do is=0,ns_ext-1
                if (is_ext_neutral(is)) cycle
                mass=mass+mp*am_ext(is)*na_ext(ix,iy,is)
                vel=vel+
     .           za_ext(ix,iy,is)*na_ext(ix,iy,is)*ua_ext(ix,iy,is)
              end do
              vth=sqrt(pz(ix,iy)/mass)
              mach(ix,iy,0)=vel/(vth*ne(ix,iy))
              mach(ix,iy,0)=pz(ix,iy)*(1.0_R8+mach(ix,iy,0)**2)
            end do
          end do
          sumjxa=0.0_R8
          sumjxi=0.0_R8
          suml=0.0_R8
          sumr=0.0_R8
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
            if (iy.ge.topcut(1)) suml(1)=suml(1)+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
            ix=nx-1
            if (nnreg(0).eq.7.and.topcut(1).lt.topcut(2)) then
              if (iy.ge.topcut(1).and.iy.lt.topcut(2))
     .              sumr(1)=sumr(1)+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
            else if (nnreg(0).ne.7) then
              if (iy.ge.topcut(1)) sumr(1)=sumr(1)+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
            end if
            if(nncut.eq.2) then
              ix=nxtr+1
              if (iy.ge.topcut(2)) suml(2)=suml(2)+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
              ix=nxtl-1
              if (nnreg(0).eq.7.and.topcut(1).lt.topcut(2)) then
                if (iy.ge.topcut(2)) sumr(2)=sumr(2)+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
              else
                if (iy.ge.topcut(1)) sumr(2)=sumr(2)+mach(ix,iy,0)*
     1                sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)*cr(ix,iy)*
     2                hy1(ix,iy)
              endif
            endif
          end do
          if (nnreg(0).eq.7.and.topcut(1).lt.topcut(2)) then
            write(*,*) 'Midplane pressure: ',sumjxa
            write(*,*) 'Inner target:      ',suml(1)
            write(*,*) 'Near SOL target:   ',sumr(1)
            write(*,*) 'Far SOL target:    ',sumr(2)
            write(*,*) 'Pressure drop (inner & outer): ',
     1             suml(1)/sumjxa,(sumr(1)+sumr(2))/sumjxa
          else if (nnreg(0).le.7) then
            write(*,*) 'Midplane pressure: ',sumjxa
            write(*,*) 'Western target:    ',suml(1)
            write(*,*) 'Eastern target:    ',sumr(1)
            write(*,*) 'Pressure drop (inner & outer): ',
     1             suml(1)/sumjxa,sumr(1)/sumjxa
          else if (nnreg(0).eq.8) then
            write(*,*) 'Outer midplane pressure: ',sumjxa
            write(*,*) 'Outer top target:        ',suml(2)
            write(*,*) 'Outer bottom target:     ',sumr(1)
            write(*,*) 'Pressure drop (top & bottom): ',
     1             suml(2)/sumjxa,sumr(1)/sumjxa
            write(*,*) ' '
            write(*,*) 'Inner midplane pressure: ',sumjxi
            write(*,*) 'Inner top target:        ',sumr(2)
            write(*,*) 'Inner bottom target:     ',suml(1)
            write(*,*) 'Pressure drop (top & bottom): ',
     1             sumr(2)/sumjxi,suml(1)/sumjxi
          endif
          call copy(mach,1.0_R8,nz,textbl)
        case ('conn')
! New treatment using connectivity arrays
          description='Connection length (m)'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              mach(ix,iy,0) =
     &         get_connection_length( 3, 2, nx, ny, ix, iy)
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('wpe')
          description='\Gx\R_{pe} (rad/s)'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              mach(ix,iy,0)=sqrt(ne(ix,iy)*qe**2/me/eps0)
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('wce')
          description='\Gx\R_{ce} (rad/s)'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              mach(ix,iy,0)=qe*bb(ix,iy,3)/me
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('wci')
          description='\Gx\R_{ci} (rad/s)'
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
          description='v_{th,e} (m/s)'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              mach(ix,iy,0)=sqrt(te(ix,iy)/me)
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('vthi')
          description='v_{th,i} (m/s)'
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
          description='v_{Alfven} (m/s)'
          nz=ns
          do is=0,ns-1
            do ix=-1,nx
              do iy=-1,ny
                if (is_neutral(is)) then
                  mach(ix,iy,is)=0.0_R8
                else
                  mach(ix,iy,is)=bb(ix,iy,3)/
     1             sqrt(mu0*na(ix,iy,is)*am(is)*mp)
                end if
              end do
            end do
          enddo
          call copy(mach,1.0_R8,nz,textpl)
        case ('cr')
          description='CR (m)'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              mach(ix,iy,0)=cr(ix,iy)
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('cz')
          description='CZ (m)'
          nz=1
          do ix=-1,nx
            do iy=-1,ny
              mach(ix,iy,0)=cz(ix,iy)
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('na0')
          description='na at previous time step (m^{-3})'
          nz=ns
          call copy(na0,1.0_R8,nz,textpl)
        case ('ne0')
          description='ne at previous time step (m^-3)'
          nz=1
          call copy(ne0,1.0_R8,nz,textbl)
        case ('ni0')
          description='ni at previous time step (m^-3)'
          nz=2
          call copy(ni0,1.0_R8,nz,textni)
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
        case ('smth')
          description='Thermal force (N)'
          nz=ns
          mach(-1:nx,-1:ny,0:nz-1)=smpr(-1:nx,-1:ny,0:nz-1)+
     .                             smpt(-1:nx,-1:ny,0:nz-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('smpr')
          description='Thermal force (proportional to exhp) (N)'
          nz=ns
          call copy(smpr,1.0_R8,nz,textpl)
        case ('smpt')
          description='Thermal force '//
     .     '(proportional to T_{i} and T_{e} gradients) (N)'
          nz=ns
          call copy(smpt,1.0_R8,nz,textpl)
        case ('smfr')
          description='Friction force (N)'
          nz=ns
          call copy(smfr,1.0_R8,nz,textpl)
        case ('smcf')
          description='Centrifugal force (N)'
          nz=ns
          call copy(smcf,1.0_R8,nz,textpl)
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
          description='Self-pressure force'
          nz=ns
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                ix1=leftix(ix,iy)
                if (ix1.ne.-2) then
                  mach(ix,iy,is)=-0.5_R8*
     1             abs(wbbl(ix,iy,0))/wbbl(ix,iy,3)*
     2             (gs(ix1,iy,0)+gs(ix,iy,0))/2.0_R8*qc(ix,iy)*
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
          stack_req=(/1,0,0,0,0,0,0/)
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
          stack_req=(/1,0,0,0,0,0,0/)
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
          stack_req=(/1,0,0,0,0,0,0/)
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
          description='Particle source (#.s^{-1})'
          nz=ns
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
            if (isel.eq.0) then
              description=  ! IGNORE
     .         'Particle source, constant part (#.s^{-1})'
            else if (isel.eq.1) then
              description=  ! IGNORE
     .         'Particle source, variable part (m^{3}.s^{-1})'
            else
              write(*,*) 'selected value of ',isel,' out of range 0..1'
              return
            end if
            do is=0,ns-1
              mach(-1:nx,-1:ny,is)=sna(-1:nx,-1:ny,isel,is)
            enddo
            write(*,*) '"sna" option with ',isel
          else
            call lower_case(description)
            description='Total ' // trim(description)  ! IGNORE
            do is=0,ns-1
              mach(-1:nx,-1:ny,is)=sna(-1:nx,-1:ny,0,is)+
     .          na(-1:nx,-1:ny,is)*sna(-1:nx,-1:ny,1,is)
            end do
          end if
          call copy(mach,1.0_R8,nz,textpl)
        case ('sne')
          description='Electron particle source (e.s^{-1})'
          nz=1
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
            if (isel.eq.0) then
              description=  ! IGNORE
     .         'Electron particle source, constant part (e.s^{-1})'
            else if (isel.eq.1) then
              description=  ! IGNORE
     .         'Electron particle source, variable part (m^{3}.s^{-1})'
            else
              write(*,*) 'selected value of ',isel,' out of range 0..1'
              return
            end if
            mach(-1:nx,-1:ny,0)=sne(-1:nx,-1:ny,isel)
            write(*,*) '"sne" option with ',isel
          else
            call lower_case(description)
            description='Total ' // trim(description)  ! IGNORE
            mach(-1:nx,-1:ny,0)=sne(-1:nx,-1:ny,0)+
     .          ne(-1:nx,-1:ny)*sne(-1:nx,-1:ny,1)
          endif
          call copy(mach,1.0_R8,nz,textbl)
        case ('smo')
          description='Parallel momentum source (N)'
          nz=ns
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
            if (isel.eq.0) then
              description=  ! IGNORE
     .         'Parallel momentum source, constant part (N)'
            else if (isel.eq.1) then
              description=  ! IGNORE
     .         'Parallel momentum source, ua part (kg.s^{-1})'
            else if (isel.eq.2) then
              description=  ! IGNORE
     .         'Parallel momentum source, ' //
     .         'mass density part (m^{4}.s^{-2})'
            else if (isel.eq.3) then
              description=  ! IGNORE
     .         'Parallel momentum source, ' //
     .         'momentum density part (m^{3}.s^{-1})'
            else
              write(*,*) 'selected value of ',isel,' out of range 0..3'
              return
            end if
            do is=0,ns-1
              mach(:,:,is)=smo(:,:,isel,is)
            enddo
            write(*,*) '"smo" option with ',isel
          else
            description='Total parallel momentum source (N)'  ! IGNORE
            do is=0,ns-1
              mach(-1:nx,-1:ny,is)=smo(-1:nx,-1:ny,0,is)+
     .          ua(-1:nx,-1:ny,is)*smo(-1:nx,-1:ny,1,is)+
     .          ro(-1:nx,-1:ny,is)*smo(-1:nx,-1:ny,2,is)+
     .          ro(-1:nx,-1:ny,is)*smo(-1:nx,-1:ny,3,is)*
     .          ua(-1:nx,-1:ny,is)
            end do
          end if
          call copy(mach,1.0_R8,nz,textpl)
        case ('smq')
          description='smq (momentum sources from atomic physics) (N)'
          nz=ns
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
            if (isel.eq.0) then
              description='smq, constant part (N)'  ! IGNORE
            else if (isel.eq.1) then
              description='smq, ua part (kg.s^{-1})'  ! IGNORE
            else if (isel.eq.2) then
              description='smq, mass density part (m^{4}.s^{-2})'  ! IGNORE
            else if (isel.eq.3) then
              description='smq, momentum density part (m^{3}.s^{-1})'  ! IGNORE
            else
              write(*,*) 'selected value of ',isel,' out of range 0..3'
              return
            end if
            do is=0,ns-1
              mach(:,:,is)=smq(:,:,isel,is)
            enddo
            write(*,*) '"smq" option with ',isel
          else
            description='Total ' // trim(description)  ! IGNORE
            do is=0,ns-1
              mach(-1:nx,-1:ny,is)=smq(-1:nx,-1:ny,0,is)+
     .          ua(-1:nx,-1:ny,is)*smq(-1:nx,-1:ny,1,is)+
     .          ro(-1:nx,-1:ny,is)*smq(-1:nx,-1:ny,2,is)+
     .          ro(-1:nx,-1:ny,is)*smq(-1:nx,-1:ny,3,is)*
     .          ua(-1:nx,-1:ny,is)
            end do
          endif
          call copy(mach,1.0_R8,nz,textpl)
        case ('she')
          description='Electron heat source (W)'
          nz=1
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
            if (isel.eq.0) then
              description=  ! IGNORE
     .         'Electron heat source, constant part (W)'
            else if (isel.eq.1) then
              description=  ! IGNORE
     .         'Electron heat source, temperature part (s^{-1})'
            else if (isel.eq.2) then
              description=  ! IGNORE
     .         'Electron heat source, density part (W.m^{3})'
            else if (isel.eq.3) then
              description=  ! IGNORE
     .         'Electron heat source, ' //
     .         'energy density part (m^{3}.s^{-1})'
            else
              write(*,*) 'selected value of ',isel,' out of range 0..3'
              return
            endif
            mach(-1:nx,-1:ny,0)=she(-1:nx,-1:ny,isel)
            write(*,*) '"she" option with ',isel
          else
            description='Total electron heat source (W)'  ! IGNORE
            mach(-1:nx,-1:ny,0)=she(-1:nx,-1:ny,0)+
     .          te(-1:nx,-1:ny)*she(-1:nx,-1:ny,1)+
     .          ne(-1:nx,-1:ny)*she(-1:nx,-1:ny,2)+
     .          ne(-1:nx,-1:ny)*she(-1:nx,-1:ny,3)*te(-1:nx,-1:ny)
          endif
          call copy(mach,1.0_R8,nz,textbl)
        case ('shi')
          description='Ion heat source (W)'
          nz=1
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
            if (isel.eq.0) then
              description=  ! IGNORE
     .         'Ion heat source, constant part (W)'
            else if (isel.eq.1) then
              description=  ! IGNORE
     .         'Ion heat source, temperature part (s^{-1})'
            else if (isel.eq.2) then
              description=  ! IGNORE
     .         'Ion heat source, density part (W.m^{3})'
            else if (isel.eq.3) then
              description=  ! IGNORE
     .         'Ion heat source, energy density part (m^{3}.s^{-1})'
            else
              write(*,*) 'selected value of ',isel,' out of range 0..3'
              return
            endif
            mach(-1:nx,-1:ny,0)=shi(-1:nx,-1:ny,isel)
            write(*,*) '"shi" option with ',isel
          else
            description='Total ion heat source (W)'  ! IGNORE
            mach(-1:nx,-1:ny,0)=shi(-1:nx,-1:ny,0)+
     .          ti(-1:nx,-1:ny)*shi(-1:nx,-1:ny,1)+
     .          ni(-1:nx,-1:ny,0)*shi(-1:nx,-1:ny,2)+
     .          ni(-1:nx,-1:ny,0)*shi(-1:nx,-1:ny,3)*ti(-1:nx,-1:ny)
          endif
          call copy(mach,1.0_R8,nz,textbl)
        case ('sch')
          description='Charge source (A)'
          nz=1
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
            if (isel.eq.0) then
              description=  ! IGNORE
     .         'Charge source, constant part (A)'
            else if (isel.eq.1) then
              description=  ! IGNORE
     .         'Charge source, potential part (S)'
            else if (isel.eq.2) then
              description=  ! IGNORE
     .         'Charge source, density part (A.m^{3})'
            else if (isel.eq.3) then
              description=  ! IGNORE
     .         'Charge source, potential*density part (S.m^{3})'
            else
              write(*,*) 'selected value of ',isel,' out of range 0..3'
              return
            endif
            mach(-1:nx,-1:ny,0)=sch(-1:nx,-1:ny,isel)
            write(*,*) '"sch" option with ',isel
          else
            description='Total charge source (A)'  ! IGNORE
            mach(-1:nx,-1:ny,0)=sch(-1:nx,-1:ny,0)+
     .          po(-1:nx,-1:ny)*sch(-1:nx,-1:ny,1)+
     .          ne(-1:nx,-1:ny)*sch(-1:nx,-1:ny,2)+
     .          ne(-1:nx,-1:ny)*sch(-1:nx,-1:ny,3)*po(-1:nx,-1:ny)
          endif
          call copy(mach,1.0_R8,nz,textbl)
        case ('smaf')
          description='smaf (extended friction term) (N)'
          nz=ns
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
            if (isel.eq.0) then
              description=  ! IGNORE
     .         'Friction term, constant part (N)'
            else if (isel.eq.1) then
              description=  ! IGNORE
     .         'Friction term, ua part (kg.s^{-1})'
            else if (isel.eq.2) then
              description=  ! IGNORE
     .         'Friction term, mass density part (m^{4}.s^{-2})'
            else if (isel.eq.3) then
              description=  ! IGNORE
     .         'Friction term, momentum density part (m^{3}.s^{-1})'
            else
              write(*,*) 'selected value of ',isel,' out of range 0..3'
              return
            endif
            do is=0,ns-1
              mach(:,:,is)=b2npmo_smaf(:,:,isel,is)
            enddo
            write(*,*)
     &         '"smaf (extended friction term)" option with ',isel
          else
            description='Total ' // trim(description)  ! IGNORE
            do is=0,ns-1
              mach(-1:nx,-1:ny,is)=b2npmo_smaf(-1:nx,-1:ny,0,is)+
     .          ua(-1:nx,-1:ny,is)*b2npmo_smaf(-1:nx,-1:ny,1,is)+
     .          ro(-1:nx,-1:ny,is)*b2npmo_smaf(-1:nx,-1:ny,2,is)+
     .          ro(-1:nx,-1:ny,is)*b2npmo_smaf(-1:nx,-1:ny,3,is)*
     .          ua(-1:nx,-1:ny,is)
            end do
          end if
          call copy(mach,1.0_R8,nz,textpl)
        case ('smag')
          description='smag (pressure gradient term) (N)'
          nz=ns
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
            if (isel.eq.0) then
              description=  ! IGNORE
     .         'Pressure gradient term, constant part (N)'
            else if (isel.eq.1) then
              description=  ! IGNORE
     .         'Pressure gradient term, ua part (kg.s^{-1})'
            else if (isel.eq.2) then
              description=  ! IGNORE
     .         'Pressure gradient term, ' //
     .         'mass density part (m^{4}.s^{-2})'
            else if (isel.eq.3) then
              description=  ! IGNORE
     .         'Pressure gradient term, ' //
     .         'momentum density part (m^{3}.s^{-1})'
            else
              write(*,*) 'selected value of ',isel,' out of range 0..3'
              return
            endif
            do is=0,ns-1
              mach(:,:,is)=b2npmo_smag(:,:,isel,is)
            enddo
            write(*,*)
     &         '"smag (pressure gradient term)" option with ',isel
          else
            description='Total ' // trim(description)  ! IGNORE
            do is=0,ns-1
              mach(-1:nx,-1:ny,is)=b2npmo_smag(-1:nx,-1:ny,0,is)+
     .          ua(-1:nx,-1:ny,is)*b2npmo_smag(-1:nx,-1:ny,1,is)+
     .          ro(-1:nx,-1:ny,is)*b2npmo_smag(-1:nx,-1:ny,2,is)+
     .          ro(-1:nx,-1:ny,is)*b2npmo_smag(-1:nx,-1:ny,3,is)*
     .          ua(-1:nx,-1:ny,is)
            end do
          end if
          call copy(mach,1.0_R8,nz,textpl)
        case ('smav')
          description='smav (additional viscosity term) (N)'
          nz=ns
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
            if (isel.eq.0) then
              description=  ! IGNORE
     .         'Additional viscosity term, constant part (N)'
            else if (isel.eq.1) then
              description=  ! IGNORE
     .         'Additional viscosity term, ua part (kg.s^{-1})'
            else if (isel.eq.2) then
              description=  ! IGNORE
     .         'Additional viscosity term, ' //
     .         'mass density part (m^{4}.s^{-2})'
            else if (isel.eq.3) then
              description=  ! IGNORE
     .         'Additional viscosity term, ' //
     .         'momentum density part (m^{3}.s^{-1})'
            else
              write(*,*) 'selected value of ',isel,' out of range 0..3'
              return
            endif
            do is=0,ns-1
              mach(:,:,is)=b2npmo_smav(:,:,isel,is)
            enddo
            write(*,*)
     &         '"smav (additional viscosity term)" option with ',isel
          else
            description='Total ' // trim(description)  ! IGNORE
            do is=0,ns-1
              mach(-1:nx,-1:ny,is)=b2npmo_smav(-1:nx,-1:ny,0,is)+
     .          ua(-1:nx,-1:ny,is)*b2npmo_smav(-1:nx,-1:ny,1,is)+
     .          ro(-1:nx,-1:ny,is)*b2npmo_smav(-1:nx,-1:ny,2,is)+
     .          ro(-1:nx,-1:ny,is)*b2npmo_smav(-1:nx,-1:ny,3,is)*
     .          ua(-1:nx,-1:ny,is)
            end do
          end if
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
          call copy(corte,1.0_R8/ev,nz,textbl)
        case ('cti')
          description='CORTI'
          nz=1
          call copy(corti,1.0_R8/ev,nz,textbl)
        case ('cpo')
          description='CORPO'
          nz=1
          call copy(corpo,1.0_R8,nz,textbl)
        case ('ixlt')
          description='ixlt'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('ixle')
          description='ixle'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(ix .le. istack(iistack)) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              end if
            end do
          end do
          write(*,*) '"ixle" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('ixeq')
          description='ixeq'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('ixne')
          description='ixne'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(ix .eq. istack(iistack)) then
                mach(ix,iy,0)=0.0_R8
              else
                mach(ix,iy,0)=1.0_R8
              end if
            end do
          end do
          write(*,*) '"ixne" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('ixge')
          description='ixge'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(ix .ge. istack(iistack)) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              end if
            end do
          end do
          write(*,*) '"ixge" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('ixgt')
          description='ixgt'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
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
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('iyle')
          description='iyle'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(iy .le. istack(iistack)) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              end if
            end do
          end do
          write(*,*) '"iyle" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('iyeq')
          description='iyeq'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('iyne')
          description='iyne'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(iy .eq. istack(iistack)) then
                mach(ix,iy,0)=0.0_R8
              else
                mach(ix,iy,0)=1.0_R8
              end if
            end do
          end do
          write(*,*) '"iyne" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('iyge')
          description='iyge'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ny
            do ix=-1,nx
              if(iy .ge. istack(iistack)) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              end if
            end do
          end do
          write(*,*) '"iyge" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,nz,textbl)
        case ('iygt')
          description='iygt'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
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
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('isle')
          description='isle'
          nz=ns
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                if(is .le. istack(iistack)) then
                  mach(ix,iy,is)=1.0_R8
                else
                  mach(ix,iy,is)=0.0_R8
                end if
              end do
            end do
          end do
          write(*,*) '"isle" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,ns,textpl)
        case ('iseq')
          description='iseq'
          nz=ns
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('isne')
          description='isne'
          nz=ns
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                if(is .eq. istack(iistack)) then
                  mach(ix,iy,is)=0.0_R8
                else
                  mach(ix,iy,is)=1.0_R8
                end if
              end do
            end do
          end do
          write(*,*) '"isne" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,ns,textpl)
        case ('isge')
          description='isge'
          nz=ns
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ns-1
            do iy=-1,ny
              do ix=-1,nx
                if(is .ge. istack(iistack)) then
                  mach(ix,iy,is)=1.0_R8
                else
                  mach(ix,iy,is)=0.0_R8
                end if
              end do
            end do
          end do
          write(*,*) '"isge" option with ',ipop(istack,iistack)
          call copy(mach,1.0_R8,ns,textpl)
        case ('isgt')
          description='isgt'
          nz=ns
          stack_req=(/1,0,0,0,0,0,0/)
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
          description='ln \GK\R'
          nz=1
          isel=0
          if(iistack.gt.0) then
            isel=ipop(istack,iistack)
          end if
          if (isel.eq.0) then
            description=trim(description)//'_{Br}'
          else if (isel.eq.1) then
            description=trim(description)//'_{ii}'
          else if (isel.eq.2) then
            description=trim(description)//'_{ee}'
          else if (isel.eq.3) then
            description=trim(description)//'_{ei}'
          else
            write(*,*) 'selected value of ',isel,' out of range 0..3'
            return
          end if
          call b2tlnl(nx,ny,isel,te,ti,ne,lnlam)
          mach(-1:nx,-1:ny,0)=lnlam(-1:nx,-1:ny)
          call copy(mach,1.0_R8,nz,textbl)
        case ('q')
          description='q'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              if(.not.on_closed_surface(ix,iy)) then
                mach(ix,iy,0)=0.0_R8
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
              if(.not.on_closed_surface(ix,iy)) then
                mach(ix,iy,0)=0.0_R8
              else
                mach(ix,iy,0)=shear(iy)
              end if
            end do
          end do
          call copy(mach,1.0_R8,nz,textbl)
        case ('alph')
          description='\Ga\R' //
     1     ' (ratio or pressure gradient to critical ballooning value)'
          nz=1
          write(*,*) 'WARNING: "alph" not yet generalised'
c     assumes that bb is OK for 1d cases
          do iy=-1,ny
            do ix=-1,nx
              if(.not.on_closed_surface(ix,iy)) then
                mach(ix,iy,0)=0.0_R8
              else
                if (bottomiy(ix,iy).ne.-2) then
                  mach(ix,iy,0)=r0*
     1             (pz(bottomix(ix,iy),bottomiy(ix,iy))-pz(ix,iy))/
     .              (rave(iy)-rave(bottomiy(ix,iy)))/
     2             (bb(ix,iy,3)**2/2.0_R8/mu0)/
     3             (0.0537_R8+0.534_R8*shear(iy)/safetyfactor(iy)**2)
                else if (topiy(ix,iy).ne.ny+1) then
                  mach(ix,iy,0)=r0*
     1             (pz(ix,iy)-pz(topix(ix,iy),topiy(ix,iy)))/
     .             (rave(topiy(ix,iy))-rave(iy))/
     2             (bb(ix,iy,3)**2/2.0_R8/mu0)/
     3             (0.0537_R8+0.534_R8*shear(iy)/safetyfactor(iy)**2)
                else
                  mach(ix,iy,0)=0.0_R8
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
          description='D_{pa} (m^{2}.s^{-1})'
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
          description='\Gv\R_{e} (m^{2}.s^{-1})'
          nz=1
          mach(:,:,0)=hce0(:,:)/ne(:,:)
          call copy(mach,1.0_R8,nz,textbl)
        case ('kyi')
          description='\Gv\R_{i} (m^{2}.s^{-1})'
          nz=ns
          mach(:,:,0:nz-1)=hcib(:,:,0:nz-1)/na(:,:,0:nz-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('kyi0')
          description='\Gv\R_{i}^{0} (m^{2}.s^{-1})'
          nz=1
          mach(:,:,0)=hci0(:,:)/ni(:,:,0)
          call copy(mach,1.0_R8,nz,textbl)
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
          description='Braams viscosity piece (m.kg^{-1}.s^{-1})'
          nz=ns
          call copy(bas_vis_piece,1.0_R8,nz,textpl)
        case ('sig')
          description='Anomalous conductivity (mho.m^{-1})'
          nz=1
          mach(:,:,0)=sig0(:,:)/(qe*ne(:,:))
          call copy(mach,1.0_R8,nz,textbl)
        case ('alf')
          description='Thermo-electric coefficient (mho.m^{-1})'
          nz=1
          mach(:,:,0)=alf0(:,:)/(ne(:,:)*sqrt(qe/te(:,:)))
          call copy(mach,1.0_R8,nz,textbl)
        case ('fchp')
          description='Parallel current (A)'
          nz=1
          call copy(fch_p(:,:,0),1.0_R8,nz,textbl)
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
          mach(:,:,0:ns-1)=uadia(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('vd')
          description='rad comp of the ' //
     *          'drift velocity [uadia.y] (m.s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=uadia(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('vadia.x')
          description='pol comp of the divergent ' //
     *     'part of the ion diamagnetic velocity [vadia.x] (m.s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=uadia(:,:,0,0:ns-1)-vaecrb(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('vadia.y')
          description='rad comp of the divergent ' //
     *     'part of the ion diamagnetic velocity [vadia.y] (m.s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=uadia(:,:,1,0:ns-1)-vaecrb(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('wadia.x')
          description='pol comp of the total ' //
     *          'ion diamagnetic velocity [wadia.x] (m.s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=wadia(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('wadia.y')
          description='rad comp of the total ' //
     *          'ion diamagnetic velocity [wadia.y] (m.s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=wadia(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('vexb.x')
          description='pol comp of the ' //
     *          'ion E\mU\dB velocity [vaecrb.x] (m.s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=vaecrb(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('vexb.y')
          description='rad comp of the ' //
     *          'ion E\mU\dB velocity [vaecrb.y] (m.s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=vaecrb(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('vedia.x')
          description='pol comp of the divergent part of ' //
     *     'the electron diamagnetic velocity [vedia.x] (m.s^{-1})'
          nz=1
          call copy(vedia(:,:,0),1.0_R8,nz,textbl)
        case ('vedia.y')
          description='rad comp of the divergent part of ' //
     *          'the electron diamagnetic velocity [vedia.y] (m.s^{-1})'
          nz=1
          call copy(vedia(:,:,1),1.0_R8,nz,textbl)
        case ('veexb.x')
          description='pol comp of the ' //
     *          'electron E\mU\dB velocity [veecrb.x] (m.s^{-1})'
          nz=1
          call copy(veecrb(:,:,0),1.0_R8,nz,textbl)
        case ('veexb.y')
          description='rad comp of the ' //
     *          'electron E\mU\dB velocity [veecrb.y] (m.s^{-1})'
          nz=1
          call copy(veecrb(:,:,1),1.0_R8,nz,textbl)
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
                  write(*,*) 'Overlapping isolating region in cell ',
     2                          ix,iy,' for isolating regions ',
     3                          nint(mach(ix,iy,0)),' and ',iso
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
            do ix = 1, bc_list_size(ibc)
              if (mach(bc_list_x(ix,ibc),bc_list_y(ix,ibc),0).ne.0) then
                write(*,*)
     1           'Overlapping boundary conditions in cell ',
     2            bc_list_x(ix,ibc),bc_list_y(ix,ibc),
     3            ' for boundary condition ',
     4            nint(mach(bc_list_x(ix,ibc),bc_list_y(ix,ibc),0)),
     5            ' and ',ibc
              else
                mach(bc_list_x(ix,ibc),bc_list_y(ix,ibc),0) =
     .           real(ibc,R8)
              end if
            end do
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('rc')
          description='Recycling conditions'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=0.0_R8
            enddo
          enddo
          do ibc=1,nstrat
            if (crcstra(ibc).eq.'S' .or. crcstra(ibc).eq.'N' .or.
     &          crcstra(ibc).eq.'W' .or. crcstra(ibc).eq.'E') then
              do ix = 1, rc_list_size(ibc)
                if (mach(rc_list_x(ix,ibc),rc_list_y(ix,ibc),0).ne.0)
     >           then
                  write(*,*)
     1             'Overlapping recycling conditions in cell ',
     2              rc_list_x(ix,ibc),rc_list_y(ix,ibc),
     3              ' for recycling condition ',
     4              nint(mach(rc_list_x(ix,ibc),rc_list_y(ix,ibc),0)),
     5              ' and ',ibc
                else
                  mach(rc_list_x(ix,ibc),rc_list_y(ix,ibc),0) =
     .             real(ibc,R8)
                end if
              end do
            end if
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('ix')
          description='ix (normalised)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(ix+1,R8)/real(nx+1,R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('iy')
          description='iy (normalised)'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(iy+1,R8)/real(ny+1,R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('iix')
          description='ix'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(ix,R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('iiy')
          description='iy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(iy,R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('lix')
          description='leftix'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(leftix(ix,iy),R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('liy')
          description='leftiy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(leftiy(ix,iy),R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('rix')
          description='rightix'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(rightix(ix,iy),R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('riy')
          description='rightiy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(rightiy(ix,iy),R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('bix')
          description='bottomix'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(bottomix(ix,iy),R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('biy')
          description='bottomiy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(bottomiy(ix,iy),R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('tix')
          description='topix'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(topix(ix,iy),R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('tiy')
          description='topiy'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0)=real(topiy(ix,iy),R8)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('rlcl')
          description='real cells'
          nz=1
          do iy=-1,ny
            do ix=-1,nx
              if(isRealCell(cflags(ix,iy,CELLFLAG_TYPE))) then
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
              if(cflags(ix,iy,CELLFLAG_TYPE) == GRID_GUARD) then
                mach(ix,iy,0)=1.0_R8
              else
                mach(ix,iy,0)=0.0_R8
              endif
            enddo
          enddo
          call copy(mach,1.0_R8,nz,textbl)
        case ('dib2')
          description='Test ion density (m^{-3})'
          if (nioni.gt.0) then
            nz=nioni
            mach(:,:,0:nz-1)=dib2(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8,nz,textin)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          endif
        case ('tib2')
          description='Test ion temperature (eV)'
          if (nioni.gt.0) then
            nz=nioni
            mach(:,:,0:nz-1)=tib2(0:nx+1,0:ny+1,1:nz,1)
            call copy(mach,1.0_R8/ev,nz,textin)
          else
            nz=1
            mach(:,:,0:nz-1)=0.0_R8
            call copy(mach,1.0_R8,nz,textbl)
          endif
        case ('brna')
          description='Recycling particle source (#.s^{-1})'
          nz=ns
          call copy(b2stbr_sna,1.0_R8,nz,textpl)
        case ('brne')
          description='Recycling electron source (e.s^{-1})'
          nz=1
          call copy(b2stbr_sne,1.0_R8,nz,textbl)
        case ('brmo')
          description='Recycling parallel momentum source (N)'
          nz=ns
          call copy(b2stbr_smo,1.0_R8,nz,textpl)
        case ('brhe')
          description='Recycling electron heat source (W)'
          nz=1
          call copy(b2stbr_she,1.0_R8,nz,textbl)
        case ('brhi')
          description='Recycling ion heat source (W)'
          nz=1
          call copy(b2stbr_shi,1.0_R8,nz,textbl)
        case ('brch')
          description='Recycling charge source (A)'
          nz=1
          call copy(b2stbr_sch,1.0_R8,nz,textbl)
        case ('bcna')
          description='Boundary particle source (#.s^{-1})'
          nz=ns
          call copy(b2stbc_sna,1.0_R8,nz,textpl)
        case ('bcne')
          description='Boundary electron source (e.s^{-1})'
          nz=1
          call copy(b2stbc_sne,1.0_R8,nz,textbl)
        case ('bcmo')
          description='Boundary parallel momentum source (N)'
          nz=ns
          call copy(b2stbc_smo,1.0_R8,nz,textpl)
        case ('bche')
          description='Boundary electron heat source (W)'
          nz=1
          call copy(b2stbc_she,1.0_R8,nz,textbl)
        case ('bchi')
          description='Boundary ion heat source (W)'
          nz=1
          call copy(b2stbc_shi,1.0_R8,nz,textbl)
        case ('bcch')
          description='Boundary charge source (A)'
          nz=1
          call copy(b2stbc_sch,1.0_R8,nz,textbl)
        case ('bmna')
          description='Additional particle source (#.s^{-1})'
          nz=ns
          call copy(b2stbm_sna,1.0_R8,nz,textpl)
        case ('bmne')
          description='Additional electron source (e.s^{-1})'
          nz=1
          call copy(b2stbm_sne,1.0_R8,nz,textbl)
        case ('bmmo')
          description='Additional parallel momentum source (N)'
          nz=ns
          call copy(b2stbm_smo,1.0_R8,nz,textpl)
        case ('bmhe')
          description='Additional electron heat source (W)'
          nz=1
          call copy(b2stbm_she,1.0_R8,nz,textbl)
        case ('bmhi')
          description='Additional ion heat source (W)'
          nz=1
          call copy(b2stbm_shi,1.0_R8,nz,textbl)
        case ('bmch')
          description='Additional charge source (A)'
          nz=1
          call copy(b2stbm_sch,1.0_R8,nz,textbl)
        case ('exna')
          description='External particle source (#.s^{-1})'
          nz=ns
          call copy(ext_sna,1.0_R8,nz,textpl)
        case ('exne')
          description='External electron source (e.s^{-1})'
          nz=1
          call copy(ext_sne,1.0_R8,nz,textbl)
        case ('exmo')
          description='External parallel momentum source (N)'
          nz=ns
          call copy(ext_smo,1.0_R8,nz,textpl)
        case ('exhe')
          description='External electron heat source (W)'
          nz=1
          call copy(ext_she,1.0_R8,nz,textbl)
        case ('exhi')
          description='External ion heat source (W)'
          nz=1
          call copy(ext_shi,1.0_R8,nz,textbl)
        case ('exch')
          description='External charge source (A)'
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
          description='Convective poloidal particle flux (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_32(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnay32')
          description='Convective radial particle flux (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_32(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnax52')
          description='Conductive poloidal particle flux (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_52(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnay52')
          description='Conductive radial particle flux (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_52(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnaxnodr')
          description='Poloidal particle flux without drifts (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_nodrift(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnaynodr')
          description='Radial particle flux without drifts (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_nodrift(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnix32')
          description='Convective poloidal total ion flux (s^{-1})'
          nz=1
          call copy(fni_32(:,:,0),1.0_R8,nz,textbl)
        case ('fniy32')
          description='Convective radial total ion flux (s^{-1})'
          nz=1
          call copy(fni_32(:,:,1),1.0_R8,nz,textbl)
        case ('fnix52')
          description='Conductive poloidal total ion flux (s^{-1})'
          nz=1
          call copy(fni_52(:,:,0),1.0_R8,nz,textbl)
        case ('fniy52')
          description='Conductive radial total ion flux (s^{-1})'
          nz=1
          call copy(fni_52(:,:,1),1.0_R8,nz,textbl)
        case ('fnex32')
          description='Convective poloidal electron flux (s^{-1})'
          nz=1
          call copy(fne_32(:,:,0),1.0_R8,nz,textbl)
        case ('fney32')
          description='Convective radial electron flux (s^{-1})'
          nz=1
          call copy(fne_32(:,:,1),1.0_R8,nz,textbl)
        case ('fnex52')
          description='Conductive poloidal electron flux (s^{-1})'
          nz=1
          call copy(fne_52(:,:,0),1.0_R8,nz,textbl)
        case ('fney52')
          description='Conductive radial electron flux (s^{-1})'
          nz=1
          call copy(fne_52(:,:,1),1.0_R8,nz,textbl)
        case ('fchx32')
          description='Convective poloidal current (A)'
          nz=1
          call copy(fch_32(:,:,0),1.0_R8,nz,textbl)
        case ('fchy32')
          description='Convective radial current (A)'
          nz=1
          call copy(fch_32(:,:,1),1.0_R8,nz,textbl)
        case ('fchx52')
          description='Conductive poloidal current (A)'
          nz=1
          call copy(fch_52(:,:,0),1.0_R8,nz,textbl)
        case ('fchy52')
          description='Conductive radial current (A)'
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
        case ('fnaxmdf')
          description='Modified poloidal particle flux (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_mdf(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnaymdf')
          description='Modified radial particle flux (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_mdf(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fhexmdf')
          description='Modified poloidal electron heat flux (W)'
          nz=1
          call copy(fhe_mdf(:,:,0),1.0_R8,nz,textbl)
        case ('fheymdf')
          description='Modified radial electron heat flux (W)'
          nz=1
          call copy(fhe_mdf(:,:,1),1.0_R8,nz,textbl)
        case ('fhixmdf')
          description='Modified poloidal ion heat flux (W)'
          nz=1
          call copy(fhi_mdf(:,:,0),1.0_R8,nz,textbl)
        case ('fhiymdf')
          description='Modified radial ion heat flux (W)'
          nz=1
          call copy(fhi_mdf(:,:,1),1.0_R8,nz,textbl)
        case ('fhexps')
          description=
     .     'Pfirsch-Schlueter poloidal electron heat flux (W)'
          nz=1
          call copy(fhePSch(:,:,0),1.0_R8,nz,textbl)
        case ('fheyps')
          description=
     .     'Pfirsch-Schlueter radial electron heat flux (W)'
          nz=1
          call copy(fhePSch(:,:,1),1.0_R8,nz,textbl)
        case ('fhixps')
          description='Pfirsch-Schlueter poloidal ion heat flux (W)'
          nz=1
          call copy(fhiPSch(:,:,0),1.0_R8,nz,textbl)
        case ('fhiyps')
          description='Pfirsch-Schlueter radial ion heat flux (W)'
          nz=1
          call copy(fhiPSch(:,:,1),1.0_R8,nz,textbl)
        case ('fnaxps')
          description=
     .     'Pfirsch-Schlueter poloidal particle flux (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fnaPSch(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnayps')
          description=
     .     'Pfirsch-Schlueter radial particle flux (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fnaPSch(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnaxfcor')
          description=
     .     'Poloidal particle flux with Coriolis term (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_fcor(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnayfcor')
          description=
     .     'Radial particle flux with Coriolis term (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_fcor(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnaxhe')
          description=
     .     'Poloidal particle flux for heat equations (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_he(:,:,0,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fnayhe')
          description=
     .     'Radial particle flux for heat equations (s^{-1})'
          nz=ns
          mach(:,:,0:ns-1)=fna_he(:,:,1,0:ns-1)
          call copy(mach,1.0_R8,nz,textpl)
        case ('fchxvisq')
          description=
     .     'Poloidal current due to heat flux viscosity tensor (A)'
          nz=1
          call copy(fchvisq(:,:,0),1.0_R8,nz,textbl)
        case ('fchyvisq')
          description=
     .     'Radial current due to heat flux viscosity tensor (A)'
          nz=1
          call copy(fchvisq(:,:,1),1.0_R8,nz,textbl)
        case ('fchxinrt')
          description=
     .     'Poloidal current due to inertia and gyroviscosity (A)'
          nz=1
          call copy(fchinert(:,:,0),1.0_R8,nz,textbl)
        case ('fchyinrt')
          description=
     .     'Radial current due to inertia and gyroviscosity (A)'
          nz=1
          call copy(fchinert(:,:,1),1.0_R8,nz,textbl)
*    neoclassical part                                               !cjk
        case ('dnneo')
          description='Neoclassical particle diffusivity D_{n} ' //
     &                '(m^{2}.s^{-1})'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(dn_neo,1.0_R8,nz,textplel)
        case ('chineo')
          description='Neoclassical heat conductivity \Gv\R ' //
     &                '(m^{2}.s^{-1})'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(chi_neo,1.0_R8,nz,textplel)
        case ('vnneo')
          description='Neoclassical particle pinch v_{n} ' //
     &                '(m.s^{-1})'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(vn_neo,1.0_R8,nz,textplel)
        case ('vwneo')
          description='Neoclassical heat pinch vw ' //
     &                '(m.s^{-1})'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(vw_neo,1.0_R8,nz,textplel)
        case ('lnn')
          description='Density gradient length (m)'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(lnn,1.0_R8,nz,textplel)
        case ('ltt')
          description='Temperature gradient length (m)'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(ltt,1.0_R8,nz,textplel)
        case ('nuab')
          description='Two species collision frequency (s^{-1})'
          nz=(ns+1)*(ns+1)
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(nuab,1.0_R8,nz,textplel2)
        case ('nua')
          description='Total collision frequency (s^{-1})'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(nua,1.0_R8,nz,textplel)
        case ('nustar')
          description='Normalized collisionality'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(nustar,1.0_R8,nz,textplel)
        case ('vth')
          description='Thermal velocity (m.s^{-1})'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(vtha,1.0_R8/1.e2,nz,textplel)
        case ('rhol')
          description='Full-B Larmor radius (m)'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(rhol,1.0_R8,nz,textplel)
        case ('deltal')
          description='Ratio of Larmor radius to minimum ' //
     &                'gradient length'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(deltal,1.0_R8,nz,textplel)
        case ('wb')
          description='Width of the banana orbit (m)'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(wb,1.0_R8,nz,textplel)
        case ('deltaw')
          description='Ratio of banana width to minimum ' //
     &                'gradient length'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(deltaw,1.0_R8,nz,textplel)
        case ('deltawn')
          description='Ratio of banana width to density ' //
     &                'gradient length'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(deltawn,1.0_R8,nz,textplel)
        case ('deltawt')
          description='Ratio of banana width to temperature ' //
     &                'gradient length'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(deltawt,1.0_R8,nz,textplel)
        case ('omega')
          description='Gyro-frequency (rad.s^{-1})'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(omega,1.0_R8,nz,textplel)
        case ('deltab')
          description='Ratio of collision frequency to ' //
     &                'gyro-frequency'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(deltab,1.0_R8,nz,textplel)
        case ('deltapar')
          description='Ratio of radial to parallel dynamics'
          nz=ns+1
          if(first_neoclassical) call initialize_neoclassical(neo_ic)
          call copy(deltapar,1.0_R8,nz,textplel)
        case ('chhichii')
          description='Chang-Hinton \Gv\R_{i} (m^{2}/s)'
          call chang_hinton(nxd,nyd,nsd_local,jsep,jxa,
     1     safetyfactor,ne,te,ti,na,mach)
          nz=ns
          call copy(mach,1.0_R8,nz,textpl)
*    end neoclassical part
        case ('crx')
          description='Corner radial coordinate (m)'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
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
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(istack(iistack).ge.0.and.istack(iistack).le.3) then
            call copy(cry(:,:,istack(iistack)),1.0_R8,nz,
     1       textcorner(istack(iistack)))
          else
            write (*,*) 'Invalid corner index value'
          endif
          call pop_stack(stack_req)
        case ('zn_ext')
          if(ns_ext.eq.0) return
          description='External species nuclear charge (e)'
          nz=ns_ext
          do iy=-1,ny
            do ix=-1,nx
              mach(ix,iy,0:ns_ext-1)=zn_ext(0:ns_ext-1)
            enddo
          enddo
          call copy(mach,1.0_R8,nz,text_ext)
        case ('za_ext')
          if(ns_ext.eq.0) return
          description='External species charge (e)'
          nz=ns_ext
          call copy(za_ext,1.0_R8,nz,text_ext)
        case ('za2_ext')
          if(ns_ext.eq.0) return
          description='Squared external species charge (e^{2})'
          nz=ns_ext
          call copy(za2_ext,1.0_R8,nz,text_ext)
        case ('na_ext')
          if(ns_ext.eq.0) return
          description='External species density (m^{-3})'
          nz=ns_ext
          call copy(na_ext,1.0_R8,nz,text_ext)
        case ('mi_ext')
          if(ns_ext.eq.0) return
          description='External species mass (kg)'
          nz=ns_ext
          do is = 0, ns_ext-1
            mach(:,:,is)=am_ext(is)
          end do
          call copy(mach,mp,nz,text_ext)
        case ('ua_ext')
          if(ns_ext.eq.0) return
          description='External species parallel velocity (m.s^{-1})'
          nz=ns_ext
          call copy(ua_ext,1.0_R8,nz,text_ext)
        case ('ta_ext')
          if(ns_ext.eq.0) return
          description='External species temperature (eV)'
          nz=ns_ext
          call copy(ta_ext,1.0_R8/ev,nz,text_ext)
        case ('fax_ext')
          if(ns_ext.eq.0) return
          description='External species poloidal flux (#.s^{-1})'
          nz=ns_ext
          mach(:,:,0:ns_ext-1)=fa_ext(:,:,0,0:ns_ext-1)
          call copy(mach,1.0_R8,nz,text_ext)
        case ('fay_ext')
          if(ns_ext.eq.0) return
          description='External species radial flux (#.s^{-1})'
          nz=ns_ext
          mach(:,:,0:ns_ext-1)=fa_ext(:,:,1,0:ns_ext-1)
          call copy(mach,1.0_R8,nz,text_ext)
        case ('ne_ext')
          description=
     .     'External species associated electron density (m^{-3})'
          nz=1
          call copy(ne_ext,1.0_R8,nz,textbl)
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
      integer i,ix,iy,is,isp,iatm

      cflag=.true.
      FVALUE: select case (cval)
        case('mdims')    ! print out matrix stack information
          do is=1,iwrk
            write(*,'(4i5,1x,a)') is,ndim(1:3,is),trim(graphlabel(is))
          enddo
        case('sumx')     ! sum over x
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call sumx()
          graphlabel(iwrk)=trim(graphlabel(iwrk))//' SUMX '
        case('sumy')     ! sum over y
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call sumy()
          graphlabel(iwrk)=trim(graphlabel(iwrk))//' SUMY '
        case('f.x')      ! line plot versus x
          stack_req=(/1,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.ny) then
            call fx(istack(iistack))
          else
            write(*,*) 'Faulty iy = ',istack(iistack),' chosen for f.x!'
          endif
          call pop_stack(stack_req)
        case('f.y')      ! line plot versus y
          stack_req=(/1,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.nx) then
            call fy(istack(iistack))
          else
            write(*,*) 'Faulty ix = ',istack(iistack),' chosen for f.y!'
          endif
          call pop_stack(stack_req)
        case('intx')     ! integrate over x
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call intx()
        case('inty')     ! integrate over y
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call inty()
        case('f.lx')     ! line plot versus log x
          stack_req=(/1,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.ny) then
            call flx(istack(iistack))
          else
            write(*,*)'Faulty iy = ',istack(iistack),' chosen for f.lx!'
          endif
          call pop_stack(stack_req)
        case('f.ly')     ! line plot versus log y
          stack_req=(/1,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.nx) then
            call fly(istack(iistack))
          else
            write(*,*)'Faulty ix = ',istack(iistack),' chosen for f.ly!'
          endif
          call pop_stack(stack_req)
        case('cont')     ! plot contours
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call cont
        case('surf')     ! plot patch plot
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call surf(1)
        case('line')     !
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call line(istack(iistack))
          call pop_stack(stack_req)
        case('psx')      ! per sx
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call psx
        case('psxperp')  ! per sxperp
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call psxperp
        case('psy')      ! per sy
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call psy
        case('abs')      ! abs()
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1                          abs(work(ixx(ix,iy),iyy(ix,iy),is,iwrk))
              end do
            end do
          end do
          graphlabel(iwrk)=trim(graphlabel(iwrk))//' (ABS)'
        case('ddx')      ! d/dx
          stack_req=(/0,0,0,1,0,0,0/)
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
          work(-1:ndim(1,iwrk),-1:ndim(2,iwrk),0:ndim(3,iwrk)-1,iwrk)=
     .     mach(-1:ndim(1,iwrk),-1:ndim(2,iwrk),0:ndim(3,iwrk)-1)
          graphlabel(iwrk)=trim(graphlabel(iwrk))//'/(m) (d/dx)'
        case('ddy')      ! d/dy
          stack_req=(/0,0,0,1,0,0,0/)
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
          work(-1:ndim(1,iwrk),-1:ndim(2,iwrk),0:ndim(3,iwrk)-1,iwrk)=
     .     mach(-1:ndim(1,iwrk),-1:ndim(2,iwrk),0:ndim(3,iwrk)-1)
          graphlabel(iwrk)=trim(graphlabel(iwrk))//'/(m) (d/dy)'
        case('sumz')     ! sum over ns
          stack_req=(/2,0,0,1,0,0,0/)
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
     1       trim(graphlabel(iwrk)) //' SUM'//
     2       trim(textwk(izstrt,iwrk))//' to'//trim(textwk(izend,iwrk))
          else
            graphlabel(iwrk)=
     1       trim(graphlabel(iwrk)) // trim(textwk(izstrt,iwrk))
          end if
          textwk(0,iwrk)=textbl(0)
          call pop_stack(stack_req)
        case('pvol')     ! per vol
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call pvol
        case('pne')      ! per ne
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          call pne
        case('m+')       ! matrix addition
          stack_req=(/0,0,0,2,0,0,0/)
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
        case('m-')       ! matrix subtraction
          stack_req=(/0,0,0,2,0,0,0/)
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
        case('m*')       ! matrix multiplication
          stack_req=(/0,0,0,2,0,0,0/)
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
        case('m/')       ! matrix division
          stack_req=(/0,0,0,2,0,0,0/)
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
        case('rm+')      ! matrix + real
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('rm-')      ! matrix - real
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('rm*')      ! matrix * real
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('rm/')      ! matrix / real
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('rm**')     ! matrix ** real
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('drop')     ! drop top of matrix stack
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return   ! 2900
          iwrk=iwrk-1
        case('dup')      ! duplicate top of matrix stack
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('swap')     ! swap top two elements of matrix stack
          stack_req=(/0,0,0,2,0,0,0/)
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
        case('over')     ! m2 m1 --> m2 m1 m2
          stack_req=(/0,0,0,2,0,0,0/)
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
        case('list')     ! list contents of matrix stack
          if(iwrk.eq.0) then
            write(*,*) 'No data on matrix stack'
            return
          end if
          do i = iwrk, 1, -1
            write(hlp_frm,'(a,i3,a)') '(',ndim(3,i)+1,'a)'
            write(*,hlp_frm) trim(graphlabel(i)),
     .       (trim(textwk(is,i)),is=0,ndim(3,i)-1)
          end do
        case('xsel')     ! select column
          stack_req=(/1,0,0,1,0,0,0/)
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
        case('ysel')     ! select row
          stack_req=(/1,0,0,1,0,0,0/)
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
        case('zsel')     ! select species
          stack_req=(/1,0,0,1,0,0,0/)
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
                arrowx(ixx(ix,iy),iyy(ix,iy),is)=
     1                      arrowx(ixx(ix,iy),iyy(ix,iy),is1)
                arrowy(ixx(ix,iy),iyy(ix,iy),is)=
     1                      arrowy(ixx(ix,iy),iyy(ix,iy),is1)
              end do
            end do
            textwk(is,iwrk)=textwk(is1,iwrk)
          end do
          ndim(3,iwrk)=1
          call pop_stack(stack_req)
        case('mmin')     ! min(m1,m2)
          stack_req=(/0,0,0,2,0,0,0/)
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
        case('mmax')     ! max(m1,m2)
          stack_req=(/0,0,0,2,0,0,0/)
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
        case('mlt')      ! m<r
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('mle')      ! m<=r
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('meq')      ! m=r
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('mne')      ! m!=r
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('mgt')      ! m>r
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('mge')      ! m>=r
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('rmmi','rmmin')    ! min(r1,m1)
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('rmma','rmmax')    ! max(r1,m1)
          stack_req=(/0,1,0,1,0,0,0/)
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
        case('rtom')     ! real --> matrix stack
          stack_req=(/0,1,0,0,0,0,0/)
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
        case('rtoe')     ! real --> Eirene stack
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
#ifdef B25_EIRENE
          if(iewrk .eq. wklng) call dropeirenestackbottom()
          iewrk=iewrk+1
          nedim(1,iewrk)=ntrii
          nedim(2,iewrk)=1
          eirene_stack(1:ntrii,0,iewrk)=rstack(irstack)
          values(38)=rstack(irstack)
          values(39)=rstack(irstack)
          textewk(0,iewrk)=textbl(0)
          write(eirene_label(iewrk),'(1pe15.8)') rpop(rstack,irstack)
#else
          write(*,*) 'Must compile b2plot with Eirene for rtoe function'
          dummy=rpop(rstack,irstack)
#endif
        case('log')      ! log10()
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('ln')       ! ln()
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('sqrt')     ! sqrt()
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('exp')      ! exp()
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('zrng')     ! select species range
          stack_req=(/2,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          if(istack(iistack-1).gt.istack(iistack)) then
            write(*,*) 'Parameters for zrng must be given in ',
     &       'ascending order!'
            return
          end if
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
                arrowx(ixx(ix,iy),iyy(ix,iy),is-istack(iistack-1))=
     1                       arrowx(ixx(ix,iy),iyy(ix,iy),is)
                arrowy(ixx(ix,iy),iyy(ix,iy),is-istack(iistack-1))=
     1                       arrowy(ixx(ix,iy),iyy(ix,iy),is)
              end do
            end do
            textwk(is-istack(iistack-1),iwrk)=textwk(is,iwrk)
          end do
          ndim(3,iwrk)=istack(iistack)-istack(iistack-1)+1
          call pop_stack(stack_req)
        case('chs')      ! change sign
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          do is=0,ndim(3,iwrk)-1
            do iy=-1,ndim(2,iwrk)
              do ix=-1,ndim(1,iwrk)
                work(ixx(ix,iy),iyy(ix,iy),is,iwrk)=
     1                       -work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              end do
            end do
          end do
          i=index(graphlabel(iwrk),' (CHS)')
          if (i.eq.0) then
           graphlabel(iwrk)=
     1      trim(graphlabel(iwrk))//' (CHS)'
          else
           write(graphlabel(iwrk),'(a,a)')
     1      graphlabel(iwrk)(1:i-1),
     2      graphlabel(iwrk)(i+6:len_trim(graphlabel(iwrk)))
          endif
        case('bala')     ! balance ???
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('traj')     ! trajectories
          stack_req=(/5,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(istack(5).gt.0) then
            stack_req=(/7+istack(iistack),0,0,0,0,0,0/)
          else
            stack_req=(/7+istack(iistack),0,0,1,0,0,0/)
          endif
          if(bad_args(stack_req)) return
          call caltraj(istack(iistack-6-istack(iistack)),
     +     7+istack(iistack))
          call pop_stack(stack_req)
        case('stoe')     ! matrix --> Eirene array
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
#ifdef B25_EIRENE
          if(iewrk .eq. wklng) call dropeirenestackbottom()
          iewrk=iewrk+1
          nedim(1,iewrk)=ntrii
          nedim(2,iewrk)=ndim(3,iwrk)
          do i = 1, ntrii
            if (b2_cell(ixtri(i),iytri(i))) then
              eirene_stack(i,0:ndim(3,iwrk)-1,iewrk)=
     .         work(ix_e2b(ixtri(i)),iytri(i)-1,0:ndim(3,iwrk)-1,iwrk)
            else
              eirene_stack(i,0:ndim(3,iwrk)-1,iewrk)=0.0_R8
            end if
          enddo
          values(38)=
     .     minval(eirene_stack(1:ntrii,0:nedim(2,iewrk)-1,iewrk))
          values(39)=
     .     maxval(eirene_stack(1:ntrii,0:nedim(2,iewrk)-1,iewrk))
          eirene_label(iewrk)=trim(graphlabel(iwrk))
          textewk(0:nedim(2,iewrk)-1,iewrk)=textwk(0:ndim(3,iwrk),iwrk)
#else
          write(*,*) 'Must compile with Eirene for stoe function'
#endif
          goto 2059
        case('stow')     ! matrix --> wall array
          stack_req=(/0,0,0,1,0,0,0/)
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
          wall_label(iwwrk)=trim(graphlabel(iwrk))
          textwwk(0:nwdim(2,iwwrk),iwwrk)=textwk(0:ndim(3,iwrk),iwrk)
          call set_wallbarlabel(wall_label(iwwrk))
          goto 2059
        case('xtow')     ! x dimension of matrix --> wall array
          stack_req=(/1,0,0,1,0,0,0/)
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
          wall_label(iwwrk)=trim(graphlabel(iwrk))
          textwwk(0:nwdim(2,iwwrk),iwwrk)=textwk(0:ndim(3,iwrk),iwrk)
          call set_wallbarlabel(wall_label(iwwrk))
          goto 2059
        case('ytow')     ! y dimension of matrix --> wall array
          stack_req=(/1,0,0,1,0,0,0/)
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
          wall_label(iwwrk)=trim(graphlabel(iwrk))
          textwwk(0:nwdim(2,iwwrk),iwwrk)=textwk(0:ndim(3,iwrk),iwrk)
          call set_wallbarlabel(wall_label(iwwrk))
          goto 2059
        case('rtow')     ! real --> wall stack
          stack_req=(/0,1,0,0,0,0,0/)
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
        case('rtot')     ! real --> target stack
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(itwrk .eq. wklng) call droptargetstackbottom()
          itwrk=itwrk+1
          target_stack(1:nbnd,1:ndepth_nml,itwrk)=rstack(irstack)
          values(32)=rstack(irstack)
          values(33)=rstack(irstack)
          write(target_label(itwrk),'(1pe15.8)') rpop(rstack,irstack)
          call set_targetbarlabel(target_label(itwrk))
        case('eirc')     ! add the Eirene atomic quantities to the B2 values
          if (eirc_done) return
          stack_req=(/0,0,0,0,0,0,0/)
          if(b2plot_debug.ge.12)
     &     write(*,*) iistack,irstack,isstack,iwrk,iwwrk,stack_req
          if(bad_args(stack_req)) write(*,*) 'Sorry!'
          do is=0,ns-1
            if (is_neutral(is)) then
              mach(:,:,is)=na(:,:,is)
              ua(:,:,is)=ua(:,:,is)*na(:,:,is)
              up(:,:,is)=up(:,:,is)*na(:,:,is)
              vv(:,:,is)=vv(:,:,is)*na(:,:,is)
              do iatm=1,natmi
                if (latmscl(iatm).ne.b2espcr(is)) cycle
                stmp2=trim(textan(iatm-1))
                if (isaletter(stmp2(len_trim(stmp2):len_trim(stmp2))))
     >           stmp2=trim(stmp2)//'0'
                write(*,'(5a,i3,a,i3,a)')
     .           'We seem to have found a match between ',
     &            trim(textpl(is)),' and ',trim(stmp2),
     &           '  (',is,' & ',iatm,')'
                do iy=-1,ny
                  do ix=-1,nx
                    mach(ix,iy,is)=mach(ix,iy,is)+dab2(ix+1,iy+1,iatm,1)
                  enddo
                enddo
                if (maxval(abs(pfluxa(:,:,iatm,1))).gt.0.0_R8) then
                  do iy=-1,ny
                    do ix=-1,nx
                      fna(ix,iy,0,is)=fna(ix,iy,0,is)+
     &                 pfluxa(ix+1,iy+1,iatm,1)*gs(ix,iy,0)*qz(ix,iy,1)
                      ua(ix,iy,is)=ua(ix,iy,is)+
     &                 pfluxa(ix+1,iy+1,iatm,1)*bb(ix,iy,0)/bb(ix,iy,3)
                      up(ix,iy,is)=up(ix,iy,is)+pfluxa(ix+1,iy+1,iatm,1)
                    enddo
                  enddo
                endif
                if (maxval(abs(rfluxa(:,:,iatm,1))).gt.0.0_R8) then
                  do iy=-1,ny
                    do ix=-1,nx
                      fna(ix,iy,1,is)=fna(ix,iy,1,is)+
     &                 rfluxa(ix+1,iy+1,iatm,1)*gs(ix,iy,1)
                      vv(ix,iy,is)=vv(ix,iy,is)+rfluxa(ix+1,iy+1,iatm,1)
                    enddo
                  enddo
                endif
              enddo
              na(:,:,is)=mach(:,:,is)
              do iy=-1,ny
                do ix=-1,nx
                  if(na(ix,iy,is).gt.0.0_R8) then
                    ua(ix,iy,is)=ua(ix,iy,is)/na(ix,iy,is)
                    up(ix,iy,is)=up(ix,iy,is)/na(ix,iy,is)
                    vv(ix,iy,is)=vv(ix,iy,is)/na(ix,iy,is)
                  else
                    ua(ix,iy,is)=0.0_R8
                    up(ix,iy,is)=0.0_R8
                    vv(ix,iy,is)=0.0_R8
                  endif
                enddo
              enddo
              fnax(:,:,is)=fna(:,:,0,is)
              fnay(:,:,is)=fna(:,:,1,is)
              arrowx(:,:,is)=up(:,:,is)
              arrowy(:,:,is)=vv(:,:,is)
              uutraj(:,:,is)=ua(:,:,is)
              vvtraj(:,:,is)=vv(:,:,is)
              ro(:,:,is)=am(is)*mp*na(:,:,is)
            endif
          enddo
          eirc_done=.true.
          return
        case('nese')     ! set ne from matrix stack
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return ! 2900
          do iy=-1,ndim(2,iwrk)
            do ix=-1,ndim(1,iwrk)
              ne(ix,iy)=work(ixx(ix,iy),iyy(ix,iy),0,iwrk)
            enddo
          enddo
          goto 2059
        case('tese')     ! set te from matrix stack
          stack_req=(/0,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          do iy=-1,ndim(2,iwrk)
            do ix=-1,ndim(1,iwrk)
              te(ix,iy)=work(ixx(ix,iy),iyy(ix,iy),0,iwrk)
            enddo
          enddo
          goto 2059        ! drop
        case('dase')     ! set dab2 from matrix stack
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('uuse')     ! set uutraj from matrix stack
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('vvse')     ! set vvtraj from matrix stack
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('arox')     ! set arrow field x-component from matrix
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('aroy')     ! set arrow field y-component from matrix
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('pos')      ! matrix with 1 if element > 0, else 0
          stack_req=(/0,0,0,1,0,0,0/)
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
        case('idrop')    ! drop top of integer stack
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp=ipop(istack,iistack)
        case('idup')     ! duplicate top of integer stack
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp=istack(iistack)
          call ipush(istack,iistack,itmp)
        case('iswap')    ! swap top two elements of integer stack
          stack_req=(/2,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call iswap(istack(iistack),istack(iistack-1))
        case('iover')    ! i2 i1 --> i2 i1 i2
          stack_req=(/2,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp=istack(iistack-1)
          call ipush(istack,iistack,itmp)
        case('rdrop')    ! drop top of real stack
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          t1=rpop(rstack,irstack)
        case('rdup')     ! duplicate top of real stack
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          t1=rstack(irstack)
          call rpush(rstack,irstack,t1)
        case('rswap')    ! swap top two elements of real stack
          stack_req=(/0,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call rswap(rstack(irstack),rstack(irstack-1))
        case('rover')    ! r2 r1 --> r2 r1 r2
          stack_req=(/0,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          t1=rstack(irstack-1)
          call rpush(rstack,irstack,t1)
        case default FVALUE
          cflag=.false.
      END SELECT FVALUE
c
      goto 2900

 2059 continue
      stack_req=(/0,0,0,1,0,0,0/)
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
      stack_req=(/0,1,0,0,0,0,0/)
      NVALUE: select case(cval)
        case('xmin')     ! Minimum x for plotting in the computational domain
          if(bad_args(stack_req)) return
          values(1)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(1)
        case('xmax')     ! Maximum x for plotting in the computational domain
          if(bad_args(stack_req)) return
          values(2)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(2)
        case('ymin')     ! Minimum y for plotting in the computational domain
          if(bad_args(stack_req)) return
          values(3)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(3)
        case('ymax')     ! Maximum y for plotting in the computational domain
          if(bad_args(stack_req)) return
          values(4)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(4)
        case('rmin')     ! Minimum r for plotting in the physical domain
          if(bad_args(stack_req)) return
          values(5)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(5)
        case('rmax')     ! Maximum r for plotting in the physical domain
          if(bad_args(stack_req)) return
          values(6)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(6)
        case('zmin')     ! Minimum z for plotting in the physical domain
          if(bad_args(stack_req)) return
          values(7)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(7)
        case('zmax')     ! Maximum z for plotting in the physical domain
          if(bad_args(stack_req)) return
          values(8)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(8)
        case('fmin')     ! Minimum function value for plotting
          if(bad_args(stack_req)) return
          values(9)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(9)
        case('fmax')     ! Maximum function value for plotting
          if(bad_args(stack_req)) return
          values(10)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(10)
        case('omin')     ! Minimum ordinate for 2-d plots
          if(bad_args(stack_req)) return
          values(11)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(11)
        case('omax')     ! Maximum ordinate for 2-d plots
          if(bad_args(stack_req)) return
          values(12)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(12)
        case('xref')     ! Reference position in x: 0.0=left, 0.5=mid, 1.0=right
          if(bad_args(stack_req)) return
          values(13)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(13)
        case('wmin')     ! wmin in pltwall
          if(bad_args(stack_req)) return
          values(14)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(14)
        case('wmax')     ! wmax in pltwall
          if(bad_args(stack_req)) return
          values(15)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(15)
        case('t_pt')     ! total flight time of particles in caltraj
          if(bad_args(stack_req)) return
          values(16)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(16)
        case('rshi')     ! relative R shift to be applied to chords
          if(bad_args(stack_req)) return
          values(17)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(17)
          rshift=values(17)
        case('zshi')     ! relative Z shift to be applied to chords
          if(bad_args(stack_req)) return
          values(18)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(18)
          zshift=values(18)
        case('r0')       ! used to specify the R position of an ellipse in chord.F
          if(bad_args(stack_req)) return
          values(19)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(19)
          r0=values(19)
        case('z0')       ! used to specify the Z position of an ellipse in chord.F
          if(bad_args(stack_req)) return
          values(20)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(20)
          z0=values(20)
        case('a0')       ! used to specify the R radius of an ellipse in chord.F
          if(bad_args(stack_req)) return
          values(21)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(21)
          a0=values(21)
        case('b0')       ! used to specify the Z radius of an ellipse in chord.F
          if(bad_args(stack_req)) return
          values(22)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(22)
          b0=values(22)
        case('tt0')      ! used to specify the value within the ellipse in chord.F
          if(bad_args(stack_req)) return
          values(23)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(23)
          tt0=values(23)
        case('trgxmin')  ! Minimum x for target plotting in the computational domain
          if(bad_args(stack_req)) return
          values(24)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(24)
        case('trgxmax')  ! Maximum x for target plotting in the computational domain
          if(bad_args(stack_req)) return
          values(25)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(25)
        case('trgymin')  ! Minimum y for target plotting in the computational domain
          if(bad_args(stack_req)) return
          values(26)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(26)
        case('trgymax')  ! Maximum y for target plotting in the computational domain
          if(bad_args(stack_req)) return
          values(27)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(27)
        case('trgrmin')  ! Minimum r for target plotting in the physical domain
          if(bad_args(stack_req)) return
          values(28)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(28)
        case('trgrmax')  ! Maximum r for target plotting in the physical domain
          if(bad_args(stack_req)) return
          values(29)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(29)
        case('trgzmin')  ! Minimum z for target plotting in the physical domain
          if(bad_args(stack_req)) return
          values(30)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(30)
        case('trgzmax')  ! Maximum z for target plotting in the physical domain
          if(bad_args(stack_req)) return
          values(31)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(31)
        case('trgfmin')  ! Minimum function value for target plotting
          if(bad_args(stack_req)) return
          values(32)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(32)
        case('trgfmax')  ! Maximum function value for target plotting
          if(bad_args(stack_req)) return
          values(33)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(33)
        case('trgomin')  ! Minimum ordinate for 2-d target plots
          if(bad_args(stack_req)) return
          values(34)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(34)
        case('trgomax')  ! Maximum ordinate for 2-d target plots
          if(bad_args(stack_req)) return
          values(35)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(35)
        case('smin')     ! Minimum ordinate for wall plots
          if(bad_args(stack_req)) return
          values(36)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(36)
        case('smax')     ! Maximum ordinate for wall plots
          if(bad_args(stack_req)) return
          values(37)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(37)
        case('emin')     ! Minimum function value for Eirene data plots
          if(bad_args(stack_req)) return
          values(38)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(38)
        case('emax')     ! Maximum function value for Eirene data plots
          if(bad_args(stack_req)) return
          values(39)=rpop(rstack,irstack)
          write(*,*) 'Set parameter ',trim(cval),' to ',values(39)
        case('neo_ic')   ! Neoclassical model number
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          neo_ic=ipop(istack,iistack)
          write(*,*) 'Set ',trim(cval),' to ', neo_ic
          rerun_neoclassical = .true.
          call initialize_neoclassical(neo_ic)
        case default NVALUE
          cflag=.false.
      END SELECT NVALUE

      return

      end subroutine select_parameters

      subroutine select_command
      integer i,in,is,ix,iy

      cflag=.true.
      stack_req=(/0,0,1,0,0,0,0/)
      SVALUE: select case(cval)
        case('extl')     ! extra label
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),stmp)
          call clear_extralabel(stmp)
        case('xlab')     ! extra label
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),stmp)
          call clear_xlabel(stmp)
        case('ylab')     ! extra label
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),stmp)
          call clear_ylabel(stmp)
        case('cltb')     ! colour table
          if(bad_args(stack_req)) return
          call set_color_map(spop(sstack,isstack))
        case('chor')     ! chord file
          if(bad_args(stack_req)) return
          call process_chords(spop(sstack,isstack))
        case('chpl')     ! plot chords
          if(bad_args(stack_req)) return
          call plot_chords(spop(sstack,isstack))
        case('cmd')      ! command file
          if(bad_args(stack_req)) return
          lu_pos=lu_pos+1
          if(lu_pos .gt. 10) then
            write(*,*) 'More than 10 nested files  --- infinite loop?'
            call xerrab ('Infinite loop command file')
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
     1           'scripts.local/commands|scripts/commands')) then
            write(*,*) 'cmd failed to open file'
            lu=save_lu(lu_pos)
            lu_pos=lu_pos-1
          else
            write(*,*) 'reading from lu ',lu
            call readline(textline,lu,ipos)
          end if
        case('cmdcdf')   ! read te, ti, ne, na, dab2, dmb2 from the b2movies cdf file and loop calling a command file for each frame
          if(bad_args(stack_req)) return
          lu_pos=lu_pos+1
          if(lu_pos .gt. 10) then
            write(*,*) 'More than 10 nested files  --- infinite loop?'
            call xerrab ('Infinite loop in cmdcdf')
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
     1           'scripts.local/commands|scripts/commands')) then
            write(*,*) 'cmd failed to open file'
            lu=save_lu(lu_pos)
            lu_pos=lu_pos-1
          else
#ifndef NO_CDF
            filename='b2movies.nc'
            call find_file(filename,ex)
            iret=nf_open(filename,NF_NOWRITE,ncidcdf)
            call check_cdf_status(iret)
            iret=nf_inq_varid(ncidcdf,'ntstep',varid)
            call check_cdf_status(iret)
            frameindex(1)=1
            iret=nf_get_var1_int(ncidcdf,varid,frameindex,ntstepcdf)
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
        case('chbx')     ! boxcar average of chords
          if(bad_args(stack_req)) return
          call process_chords_boxcar(spop(sstack,isstack))
        case('chvl')     ! value at positions
          if(bad_args(stack_req)) return
          call process_chord_slice(spop(sstack,isstack))
        case('flux')     ! call WOSlib flux routine
          if(bad_args(stack_req)) return
          call flux(spop(sstack,isstack),nx,ny,
     1           values(5),values(6),values(7),values(8))
        case('chmp')     ! plots a contour plot of the chord contribution to the mesh element
          if(bad_args(stack_req)) return
          call map_chords(spop(sstack,isstack))
        case('cam')      ! camera plot
          if(bad_args(stack_req)) return
          call process_camera(spop(sstack,isstack))
        case('movi')     ! movie
          if(bad_args(stack_req)) return
          moviename=spop(sstack,isstack)
        case('sonn')     ! something to do with Sonnet (no longer working)
          if(bad_args(stack_req)) return
c          sonnetfile=spop(sstack,isstack)
c          write(*,*) 'Reading a new geometry'
c          iuserinit=0
c          read_xdr_geometry=.true.
c          call b2igeo
c          write(*,*) 'Read a new geometry'
          stmp=spop(sstack,isstack)
          write(*,*) 'Read a new Sonnet geometry : No longer supported!'
        case('saw')      ! reads external data onto the matrix stack
          stack_req=(/2,0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          arg=spop(sstack,isstack)
          ny=ipop(istack,iistack)-2
          nx=ipop(istack,iistack)-2
          nz=1
          open(99,file=arg,status='old')
          do ix=-1,nx
            do iy=-1,ny
              read(99,*) dummy,dummy,
     1                   crx(ix,iy,1),cry(ix,iy,1),
     2                   crx(ix,iy,3),cry(ix,iy,3),
     3                   crx(ix,iy,0),cry(ix,iy,0),
     4                   crx(ix,iy,2),cry(ix,iy,2),
     5                   mach(ix,iy,0)
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
          values(1)=minval(xx)                                   !xmin
          values(2)=maxval(xx)+0.5_R8                            !xmax
          values(3)=minval(yy)                                   !ymin
          values(4)=maxval(yy)+0.5_R8                            !ymax
          values(5)=minval(crx)                                  !rmin
          values(6)=maxval(crx)                                  !rmax
          values(7)=minval(cry)                                  !zmin
          values(8)=maxval(cry)                                  !zmax
          write(*,*) 'values',(i,values(i),i=1,8)
          call copy(mach,1.0_R8,nz,textbl)
          graphlabel(iwrk)=arg
        case('saws')     ! reads external data about plasma species onto the matrix stack
          stack_req=(/3,0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          arg=spop(sstack,isstack)
          nz=ipop(istack,iistack)
          ny=ipop(istack,iistack)-2
          nx=ipop(istack,iistack)-2
          open(99,file=arg,status='old')
          do ix=-1,nx
            do iy=-1,ny
              read(99,*) dummy,dummy,
     1                   crx(ix,iy,1),cry(ix,iy,1),
     2                   crx(ix,iy,3),cry(ix,iy,3),
     3                   crx(ix,iy,0),cry(ix,iy,0),
     4                   crx(ix,iy,2),cry(ix,iy,2),
     5                  (mach(ix,iy,is),is=0,nz-1)
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
          values(1)=minval(xx)                                   !xmin
          values(2)=maxval(xx)+0.5_R8                            !xmax
          values(3)=minval(yy)                                   !ymin
          values(4)=maxval(yy)+0.5_R8                            !ymax
          values(5)=minval(crx)                                  !rmin
          values(6)=maxval(crx)                                  !rmax
          values(7)=minval(cry)                                  !zmin
          values(8)=maxval(cry)                                  !zmax
          write(*,*) 'values',(i,values(i),i=1,8)
          call copy(mach,1.0_R8,nz,textpl)
          graphlabel(iwrk)=arg
        case('sawx')     ! reads external data about external species onto the matrix stack
          stack_req=(/3,0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          arg=spop(sstack,isstack)
          nz=ipop(istack,iistack)
          ny=ipop(istack,iistack)-2
          nx=ipop(istack,iistack)-2
          open(99,file=arg,status='old')
          do ix=-1,nx
            do iy=-1,ny
              read(99,*) dummy,dummy,
     1                   crx(ix,iy,1),cry(ix,iy,1),
     2                   crx(ix,iy,3),cry(ix,iy,3),
     3                   crx(ix,iy,0),cry(ix,iy,0),
     4                   crx(ix,iy,2),cry(ix,iy,2),
     5                  (mach(ix,iy,is),is=0,nz-1)
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
          values(1)=minval(xx)                                   !xmin
          values(2)=maxval(xx)+0.5_R8                            !xmax
          values(3)=minval(yy)                                   !ymin
          values(4)=maxval(yy)+0.5_R8                            !ymax
          values(5)=minval(crx)                                  !rmin
          values(6)=maxval(crx)                                  !rmax
          values(7)=minval(cry)                                  !zmin
          values(8)=maxval(cry)                                  !zmax
          write(*,*) 'values',(i,values(i),i=1,8)
          call copy(mach,1.0_R8,nz,text_ext)
          graphlabel(iwrk)=arg
        case('oson')     ! output geometry in Sonnet format
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
        case('uedgegrd') ! output geometry in UEDGE format (string + real arg)
          stack_req=(/0,1,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          arg=spop(sstack,isstack)
          t1=rpop(rstack,irstack)
          open(99,file=arg)
cdpc still in testing phase
          write(99,'(5i4)') nx,ny,leftcut(1),rightcut(1),topcut(1)
          write(99,'()')         ! rm
          write(99,'(1p,3e23.15)') ((cr(ix,iy),ix=-1,nx),iy=-1,ny),
     1     ((crx(ix,iy,0),ix=-1,nx),iy=-1,ny),
     1     ((crx(ix,iy,1),ix=-1,nx),iy=-1,ny),
     1     ((crx(ix,iy,2),ix=-1,nx),iy=-1,ny),
     1     ((crx(ix,iy,3),ix=-1,nx),iy=-1,ny)
          write(99,'()')         ! zm
          write(99,'(1p,3e23.15)') ((cz(ix,iy)+t1,ix=-1,nx),iy=-1,ny),
     1     ((cry(ix,iy,0)+t1,ix=-1,nx),iy=-1,ny),
     1     ((cry(ix,iy,1)+t1,ix=-1,nx),iy=-1,ny),
     1     ((cry(ix,iy,2)+t1,ix=-1,nx),iy=-1,ny),
     1     ((cry(ix,iy,3)+t1,ix=-1,nx),iy=-1,ny)
          write(99,'()')         ! psi
          uedge_c=0.0_R8
          uedge_v=0.0_R8
          write(99,'(1p,3e23.15)') ((uedge_c(ix,iy),ix=-1,nx),iy=-1,ny),
     1     (((uedge_v(ix,iy,in),ix=-1,nx),iy=-1,ny),in=1,4)
          write(99,'()')         ! br
          uedge_c=0.0_R8
          uedge_v=0.0_R8
          write(99,'(1p,3e23.15)') ((uedge_c(ix,iy),ix=-1,nx),iy=-1,ny),
     1     (((uedge_v(ix,iy,in),ix=-1,nx),iy=-1,ny),in=1,4)
          write(99,'()')         ! bz
          uedge_c=0.0_R8
          uedge_v=0.0_R8
          write(99,'(1p,3e23.15)') ((uedge_c(ix,iy),ix=-1,nx),iy=-1,ny),
     1     (((uedge_v(ix,iy,in),ix=-1,nx),iy=-1,ny),in=1,4)
          write(99,'()')         ! bpol
          uedge_c(:,:)=bb(:,:,0)
          write(*,*) 'bpol_c: min & max ',
     1     minval(uedge_c),maxval(uedge_c)
          uedge_v=0.0_R8
          call my_intcor(nx,ny,vol,periodic_bc,nncut,
     1     leftix,leftiy,rightix,rightiy,
     1     bottomix,bottomiy,topix,topiy,region,uedge_c,corners)
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
          write(99,'()')         ! bphi
          uedge_c(:,:)=bb(:,:,2)
          write(*,*) 'bphi_c: min & max ',
     1     minval(uedge_c),maxval(uedge_c)
          uedge_v=0.0_R8
          call my_intcor(nx,ny,vol,periodic_bc,nncut,
     1     leftix,leftiy,rightix,rightiy,
     1     bottomix,bottomiy,topix,topiy,region,uedge_c,corners)
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
          write(99,'()')         ! b
          uedge_c(:,:)=bb(:,:,3)
          write(*,*) 'b_c: min & max ',
     1     minval(uedge_c),maxval(uedge_c)
          uedge_v=0.0_R8
          call my_intcor(nx,ny,vol,periodic_bc,nncut,
     1     leftix,leftiy,rightix,rightiy,
     1     bottomix,bottomiy,topix,topiy,region,uedge_c,corners)
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
        case('expf')     ! add experimental data to plot
          if(bad_args(stack_req)) return
          expfile=spop(sstack,isstack)
          call plot_expfile
        case('outp')     ! prints the value of a specific element
          stack_req=(/2,0,1,1,0,0,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          iy=ipop(istack,iistack)
          if(-1.gt.iy.or.iy.gt.ndim(2,iwrk)) then
            write(*,*) 'Faulty parameter iy = ',iy,' for outp!'
            return
          endif
          ix=ipop(istack,iistack)
          if(-1.gt.ix.or.ix.gt.ndim(1,iwrk)) then
            write(*,*) 'Faulty parameter ix = ',ix,' for outp!'
            return
          endif
          do is=0,ndim(3,iwrk)-1
            if(ndim(3,iwrk) .eq. 1) then
              write(*,'(a,1x,1p,g15.8)')
     1                 trim(stmp)//trim(textwk(is,iwrk)),
     1                 work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              if(writing_to_file) then
                write(nwtf,'(a,1x,1p,g15.8)')
     1                 trim(stmp)//trim(textwk(is,iwrk)),
     1                 work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              endif
            else
              write(*,'(a,i3,1x,1p,g15.8)')
     1                 trim(stmp)//trim(textwk(is,iwrk)),is,
     1                 work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              if(writing_to_file) then
                write(nwtf,'(a,i3,1x,1p,g15.8)')
     1                 trim(stmp)//trim(textwk(is,iwrk)),is,
     1                 work(ixx(ix,iy),iyy(ix,iy),is,iwrk)
              endif
            endif
          enddo
        case('pdev')     ! PLPLOT device name
          if(bad_args(stack_req)) return
          call set_device(spop(sstack,isstack))
        case('mprt')     ! output the top of matrix stack to a file
          stack_req=(/0,0,1,1,0,0,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          open(99,file=stmp)
          do is=0,ndim(3,iwrk)-1
            write(99,'(''#'',1x,a,1x,a)')
     1       trim(graphlabel(iwrk)),trim(textwk(is,iwrk))
            do iy=-1,ny
              do ix=-1,nx
                write(99,*) ix,iy,is,work(ix,iy,is,iwrk)
              enddo
            enddo
          enddo
          close(99)
        case('mtab')     ! output a table from the top of matrix stack to a file
          stack_req=(/4,0,1,1,0,0,0/)
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
          write(hlp_frm,'(a,i3,a)') '(a1,1x,a,1x,',nz,'a13)'
          write(99,hlp_frm) '#',
     1       trim(graphlabel(iwrk)),(trim(textwk(is,iwrk)),is=0,nz-1)
          write(hlp_frm,'(a,i3,a)') '(2i4,1p,',nz,'(1x,1e12.5))'
          do iy=iy1,iy2
            do ix=ix1,ix2
              write(99,hlp_frm) ix,iy,(work(ix,iy,i,iwrk),i=0,nz-1)
            enddo
          enddo
          close(99)
        case('amju','amjuel')    ! accesses Eirene type atomic physics data
          stack_req=(/0,0,4,0,0,0,0/)
          if(bad_args(stack_req)) return
          stmp = spop(sstack,isstack)
          call strip_spaces(stmp)
          call strcopy(3,256,crc,stmp)
          stmp = spop(sstack,isstack)
          call strip_spaces(stmp)
          call strcopy(9,256,reac,stmp)
          stmp = spop(sstack,isstack)
          call strip_spaces(stmp)
          call strcopy(4,256,h123,stmp)
          stmp = spop(sstack,isstack)
          call strip_spaces(stmp)
          call strcopy(8,256,filnam,stmp)
          nz=1
          call get_amjuel_data(filnam,h123,reac,crc,
     1     nxd,nyd,nx,ny,ne,te,mach)
          call copy(mach,1.0_R8,1,textbl)
          ndim(3,iwrk)=nz
          graphlabel(iwrk)=
     1     trim(filnam)//' '//trim(h123)//' '//trim(reac)//' '//crc
        case('adas')     ! accesses ADAS atomic physics data
          stack_req=(/0,0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          stmp = spop(sstack,isstack)
          call strip_spaces(stmp)
          call strcopy(9,256,reac,stmp)
          nz=ns
          call get_adas_data(reac,nz,nxd,nyd,nx,ny,ne,te,mach)
          call copy(mach,1.0_R8,nz,textpl)
          ndim(3,iwrk)=nz
          graphlabel(iwrk)='ADAS '//trim(reac)//' rate'
        case('scat')     ! concatenate two strings on top of string stack
          stack_req=(/0,0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          stmp=trim(spop(sstack,isstack)) // stmp
          call spush(sstack,isstack,stmp)
        case('sdrop')    ! drop top of string stack
          stmp=spop(sstack,isstack)
        case('sdup')     ! duplicate top of string stack
          stmp=sstack(isstack)
          call spush(sstack,isstack,stmp)
        case('sswap')    ! swap top two elements of string stack
          stack_req=(/0,0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          stmp1=spop(sstack,isstack)
          stmp2=spop(sstack,isstack)
          call spush(sstack,isstack,stmp1)
          call spush(sstack,isstack,stmp2)
        case('sover')    ! str2 str1 --> str2 str1 str2
          stack_req=(/0,0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          stmp=sstack(isstack-1)
          call spush(sstack,isstack,stmp)
        case('sprt')     ! output the string stack
          stack_req=(/0,0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is = isstack, 1, -1
            write(*,'(a)') trim(sstack(is))
            if(writing_to_file) then
              write(nwtf,'(a)') trim(sstack(is))
            endif
          enddo
        case('omovie')   ! old movie option
          stack_req=(/2,0,2,0,0,0,0/)
          if(bad_args(stack_req)) return
          varname=spop(sstack,isstack)
          filename=spop(sstack,isstack)
          iend=ipop(istack,iistack)
          istart=ipop(istack,iistack)
          itmp=0
#ifndef NO_CDF
          call get_cdf_ntstep(filename,itmp,ncidcdf)
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
            call get_cdf(filename, ncidcdf, is, movietime(1),
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
          stack_req=(/2,0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          varname=spop(sstack,isstack)
          filename='b2movies.nc'
          call find_file(filename,ex)
          iend=ipop(istack,iistack)
          istart=ipop(istack,iistack)
          itmp=0
#ifndef NO_CDF
          call get_cdf_ntstep(filename,itmp,ncidcdf)
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
            call get_cdf(filename, ncidcdf, is, movietime(1),
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
        case('wallblbl') ! set wall bar label
          if(bad_args(stack_req)) return
          call set_wallbarlabel(spop(sstack,isstack))
        case('trgblbl')  ! set target bar label
          if(bad_args(stack_req)) return
          call set_targetbarlabel(spop(sstack,isstack))
        case('glbl')     ! set graph label
          stack_req=(/0,0,1,1,0,0,0/)
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),graphlabel(iwrk))
        case('elbl')     ! set Eirene graph label
          stack_req=(/0,0,1,0,1,0,0/)
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),eirene_label(iewrk))
        case('wlbl')     ! set wall graph label
          stack_req=(/0,0,1,0,0,1,0/)
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),wall_label(iwwrk))
        case('tlbl')     ! set target graph label
          stack_req=(/0,0,1,0,0,0,1/)
          if(bad_args(stack_req)) return
          call customise_label(spop(sstack,isstack),target_label(itwrk))
        case('slbl')     ! set species label
          stack_req=(/1,0,1,1,0,0,0/)
          if(bad_args(stack_req)) return
          is1=ipop(istack,iistack)
          if (0.gt.is1.or.is1.ge.ndim(3,iwrk)) then
            write(*,*) 'Faulty parameter is = ',is1,' for slbl!'
            return
          endif
          call strcopy(16,256,textwk(is1,iwrk),spop(sstack,isstack))
        case('eslbl')    ! set Eirene species label
          stack_req=(/1,0,1,0,1,0,0/)
          if(bad_args(stack_req)) return
          is1=ipop(istack,iistack)
          if (0.gt.is1.or.is1.ge.nedim(2,iewrk)) then
            write(*,*) 'Faulty parameter is = ',is1,' for eslbl!'
            return
          endif
          call strcopy(16,256,textewk(is1,iewrk),spop(sstack,isstack))
        case('wslbl')    ! set wall species label
          stack_req=(/1,0,1,0,0,1,0/)
          if(bad_args(stack_req)) return
          is1=ipop(istack,iistack)
          if (0.gt.is1.or.is1.ge.nwdim(2,iwwrk)) then
            write(*,*) 'Faulty parameter is = ',is1,' for wslbl!'
            return
          endif
          call strcopy(16,256,textwwk(is1,iwwrk),spop(sstack,isstack))
        case('tranprof') ! output transport profile
          stack_req=(/0,0,0,0,0,0,0/)
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

      subroutine select_eirene_stack

#ifdef B25_EIRENE
      integer i, is

      cflag=.true.
      EVALUE: select case(cval)
        case('pdena')
          description='Atoms particle density (m^{-3})'
          if (natmi.gt.0) then
            nz=natmi
            call eirene_copy(pdena,1.0e6_R8,nz,textan)
          else
            nz=1
            call eirene_copy(pdena,1.0_R8,nz,textbl)
          end if
        case('pdenm')
          description='Molecules particle density (m^{-3})'
          if (nmoli.gt.0) then
            nz=nmoli
            call eirene_copy(pdenm,1.0e6_R8,nz,textmn)
          else
            nz=1
            call eirene_copy(pdenm,1.0_R8,nz,textbl)
          end if
        case('pdeni')
          description='Test ions particle density (m^{-3})'
          if (nioni.gt.0) then
            nz=nioni
            call eirene_copy(pdeni,1.0e6_R8,nz,textin)
          else
            nz=1
            call eirene_copy(pdeni,1.0_R8,nz,textbl)
          end if
        case('edena')
          description='Atoms energy density (eV.m^{-3})'
          if (natmi.gt.0) then
            nz=natmi
            call eirene_copy(edena,1.0e6_R8,nz,textan)
          else
            nz=1
            call eirene_copy(edena,1.0_R8,nz,textbl)
          endif
        case('edenm')
          description='Molecules energy density (eV.m^{-3})'
          if (nmoli.gt.0) then
            nz=nmoli
            call eirene_copy(edenm,1.0e6_R8,nz,textmn)
          else
            nz=1
            call eirene_copy(edenm,1.0_R8,nz,textbl)
          end if
        case('edeni')
          description='Test ions energy density (eV.m^{-3})'
          if (nioni.gt.0) then
            nz=nioni
            call eirene_copy(edeni,1.0e6_R8,nz,textin)
          else
            nz=1
            call eirene_copy(edeni,1.0_R8,nz,textbl)
          end if
        case('massa')
          description='Atom mass (kg)'
          if (natmi.gt.0) then
            nz=natmi
            do is = 1, nz
              do i = 1, ntrii
                emach(is,i) = nmassa(is)*mp
              end do
            end do
            call eirene_copy(emach(1:nz,1:ntrii),1.0_R8,nz,textan)
          else
            nz=1
            emach(1,1:ntrii) = 0.0_R8
            call eirene_copy(emach(1:nz,1:ntrii),1.0_R8,nz,textbl)
          end if
        case('massm')
          description='Molecule mass (kg)'
          if (nmoli.gt.0) then
            nz=nmoli
            do is = 1, nz
              do i = 1, ntrii
                emach(is,i) = nmassm(is)*mp
              end do
            end do
            call eirene_copy(emach(1:nz,1:ntrii),1.0_R8,nz,textmn)
          else
            nz=1
            emach(1,1:ntrii) = 0.0_R8
            call eirene_copy(emach(1:nz,1:ntrii),1.0_R8,nz,textbl)
          end if
        case('massi')
          description='Molecular ion mass (kg)'
          if (nioni.gt.0) then
            nz=nioni
            do is = 1, nz
              do i = 1, ntrii
                emach(is,i) = nmassi(is)*mp
              end do
            end do
            call eirene_copy(emach(1:nz,1:ntrii),1.0_R8,nz,textin)
          else
            nz=1
            emach(1,1:ntrii) = 0.0_R8
            call eirene_copy(emach(1:nz,1:ntrii),1.0_R8,nz,textbl)
          end if
        case('vxa')
          description=
     .     'Atoms Cartesian X-directed velocity (m.s^{-1})'
          if (natmi.gt.0) then
            nz=natmi
            call eirene_copy(vxa,1.0_R8,nz,textan)
          else
            nz=1
            call eirene_copy(vxa,1.0_R8,nz,textbl)
          end if
        case('vxm')
          description=
     .     'Molecules Cartesian X-directed velocity (m.s^{-1})'
          if (nmoli.gt.0) then
            nz=nmoli
            call eirene_copy(vxm,1.0_R8,nz,textmn)
          else
            nz=1
            call eirene_copy(vxm,1.0_R8,nz,textbl)
          end if
        case('vxi')
          description=
     .     'Test ions Cartesian X-directed velocity (m.s^{-1})'
          if (nioni.gt.0) then
            nz=nioni
            call eirene_copy(vxi,1.0_R8,nz,textin)
          else
            nz=1
            call eirene_copy(vxi,1.0_R8,nz,textbl)
          end if
        case('vya')
          description=
     .     'Atoms Cartesian Y-directed velocity (m.s^{-1})'
          if (natmi.gt.0) then
            nz=natmi
            call eirene_copy(vya,1.0_R8,nz,textan)
          else
            nz=1
            call eirene_copy(vya,1.0_R8,nz,textbl)
          end if
        case('vym')
          description=
     .     'Molecules Cartesian Y-directed velocity (m.s^{-1})'
          if (nmoli.gt.0) then
            nz=nmoli
            call eirene_copy(vym,1.0_R8,nz,textmn)
          else
            nz=1
            call eirene_copy(vym,1.0_R8,nz,textbl)
          end if
        case('vyi')
          description=
     .     'Test ions Cartesian Y-directed velocity (m.s^{-1})'
          if (nioni.gt.0) then
            nz=nioni
            call eirene_copy(vyi,1.0_R8,nz,textin)
          else
            nz=1
            call eirene_copy(vyi,1.0_R8,nz,textbl)
          end if
        case('vza')
          description=
     .     'Atoms Cartesian Z-directed velocity (m.s^{-1})'
          if (natmi.gt.0) then
            nz=natmi
            call eirene_copy(vza,1.0_R8,nz,textan)
          else
            nz=1
            call eirene_copy(vza,1.0_R8,nz,textbl)
          end if
        case('vzm')
          description=
     .     'Molecules Cartesian Z-directed velocity (m.s^{-1})'
          if (nmoli.gt.0) then
            nz=nmoli
            call eirene_copy(vzm,1.0_R8,nz,textmn)
          else
            nz=1
            call eirene_copy(vzm,1.0_R8,nz,textbl)
          end if
        case('vzi')
          description=
     .     'Test ions Cartesian Z-directed velocity (m.s^{-1})'
          if (nioni.gt.0) then
            nz=nioni
            call eirene_copy(vzi,1.0_R8,nz,textin)
          else
            nz=1
            call eirene_copy(vzi,1.0_R8,nz,textbl)
          end if
        case('vxdena')
          description=
     .     'Atoms Cartesian X-momentum density (kg.m^{-2}.s^{-1})'
          if (natmi.gt.0) then
            nz=natmi
            call eirene_copy(vxdena,10.0_R8,nz,textan)
          else
            nz=1
            call eirene_copy(vxdena,1.0_R8,nz,textbl)
          end if
        case('vxdenm')
          description=
     .     'Molecules Cartesian X-momentum density (kg.m^{-2}.s^{-1})'
          if (nmoli.gt.0) then
            nz=nmoli
            call eirene_copy(vxdenm,10.0_R8,nz,textmn)
          else
            nz=1
            call eirene_copy(vxdenm,1.0_R8,nz,textbl)
          end if
        case('vxdeni')
          description=
     .     'Test ions Cartesian X-momentum density (kg.m^{-2}.s^{-1})'
          if (nioni.gt.0) then
            nz=nioni
            call eirene_copy(vxdeni,10.0_R8,nz,textin)
          else
            nz=1
            call eirene_copy(vxdeni,1.0_R8,nz,textbl)
          end if
        case('vydena')
          description=
     .     'Atoms Cartesian Y-momentum density (kg.m^{-2}.s^{-1})'
          if (natmi.gt.0) then
            nz=natmi
            call eirene_copy(vydena,10.0_R8,nz,textan)
          else
            nz=1
            call eirene_copy(vydena,1.0_R8,nz,textbl)
          end if
        case('vydenm')
          description=
     .     'Molecules Cartesian Y-momentum density (kg.m^{-2}.s^{-1})'
          if (nmoli.gt.0) then
            nz=nmoli
            call eirene_copy(vydenm,10.0_R8,nz,textmn)
          else
            nz=1
            call eirene_copy(vydenm,1.0_R8,nz,textbl)
          end if
        case('vydeni')
          description=
     .     'Test ions Cartesian Y-momentum density (kg.m^{-2}.s^{-1})'
          if (nioni.gt.0) then
            nz=nioni
            call eirene_copy(vydeni,10.0_R8,nz,textin)
          else
            nz=1
            call eirene_copy(vydeni,1.0_R8,nz,textbl)
          end if
        case('vzdena')
          description=
     .     'Atoms Cartesian Z-momentum density (kg.m^{-2}.s^{-1})'
          if (natmi.gt.0) then
            nz=natmi
            call eirene_copy(vzdena,10.0_R8,nz,textan)
          else
            nz=1
            call eirene_copy(vzdena,1.0_R8,nz,textbl)
          end if
        case('vzdenm')
          description=
     .     'Molecules Cartesian Z-momentum density (kg.m^{-2}.s^{-1})'
          if (nmoli.gt.0) then
            nz=nmoli
            call eirene_copy(vzdenm,10.0_R8,nz,textmn)
          else
            nz=1
            call eirene_copy(vzdenm,1.0_R8,nz,textbl)
          end if
        case('vzdeni')
          description=
     .     'Test ions Cartesian Z-momentum density (kg.m^{-2}.s^{-1})'
          if (nioni.gt.0) then
            nz=nioni
            call eirene_copy(vzdeni,10.0_R8,nz,textin)
          else
            nz=1
            call eirene_copy(vzdeni,1.0_R8,nz,textbl)
          end if
        case('triarea')
          description='Triangle area (m^{2})'
          nz=1
          call eirene_copy(triangle_area,1.0_R8,nz,textbl)
        case('triavol')
          description='Triangle volume (m^{3})'
          nz=1
          call eirene_copy(triangle_vol,1.0_R8,nz,textbl)
        case('tr')
          description='Triangle centroid R coordinate (m)'
          nz=1
          do i = 1, ntrii
            emach(1,i) = (xtrian(necke(1,i))+xtrian(necke(2,i))+
     .                    xtrian(necke(3,i)))/3.0_R8
          end do
          call eirene_copy(emach(1:nz,1:ntrii),1.0_R8,nz,textbl)
        case('tz')
          description='Triangle centroid Z coordinate (m)'
          nz=1
          do i = 1, ntrii
            emach(1,i) = (ytrian(necke(1,i))+ytrian(necke(2,i))+
     .                    ytrian(necke(3,i)))/3.0_R8
          end do
          call eirene_copy(emach(1:nz,1:ntrii),1.0_R8,nz,textbl)
        case('pux')
          description='PUX'
          nz=1
          call eirene_coarse_copy(pux,1.0_R8,nz,textbl)
        case('puy')
          description='PUY'
          nz=1
          call eirene_coarse_copy(puy,1.0_R8,nz,textbl)
        case('pvx')
          description='PVX'
          nz=1
          call eirene_coarse_copy(pvx,1.0_R8,nz,textbl)
        case('pvy')
          description='PVY'
          nz=1
          call eirene_coarse_copy(pvy,1.0_R8,nz,textbl)
        case('tindex')
          description='Triangle index'
          nz=1
          do i = 1, ntrii
            emach(1,i) = real(i)
          end do
          call eirene_copy(emach(1:nz,1:ntrii),1.0_R8,nz,textbl)
        case default EVALUE
          cflag=.false.
      END SELECT EVALUE

      if(cflag) then
        call customise_label(description,eirene_label(iewrk))
        call set_wallbarlabel(eirene_label(iewrk))
        if (b2plot_debug.ge.1)
     .   write(*,*) 'Copied ',trim(cval),' to Eirene stack.'
      endif
      return
#endif
      end subroutine select_eirene_stack

      subroutine select_eirene_command
#ifdef B25_EIRENE
      integer i, ix, iy, is, it, iclv
      integer :: nntria(0:nxd-1,0:nyd-1)
      real (kind=R8) :: volsum,
     *   xc, yc, x1, x2, x3, x4, y1, y2, y3, y4
      real (kind=R8), allocatable :: eavg(:)
      logical inside_rectangle
      external inside_rectangle

      cflag=.true.
      EFVALUE: select case (cval)
        case('edims')    ! print out Eirene stack information
          do is=1,iewrk
            write(*,'(3i5,1x,a)')
     .       is,nedim(1:2,is),trim(eirene_label(is))
          enddo
        case('psta')     ! per triangle area
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call psta
        case('pstv')     ! per triangle volume
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call pstv
        case('sumt')     ! sum over triangles
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call sumt()
          eirene_label(iewrk)=trim(eirene_label(iewrk))//' SUM '
        case('esurf')    ! plot triangle patch plot
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          call triangle_surf()
        case('eclvs')    ! customise Eirene contour levels
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          stack_req=(/1,istack(iistack),0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call eclvs(istack(iistack),rstack,irstack)
          do iclv=1,istack(iistack)
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('eclvd')    ! customise Eirene contour level delta
          stack_req=(/0,3,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call eclvd(rstack(irstack-2),rstack(irstack-1),
     1               rstack(irstack))
          do iclv=1,3
            dummy=rpop(rstack,irstack)
          end do
        case('eclvn')    ! customise Eirene contour level range and number
          stack_req=(/1,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call eclvn(rstack(irstack-1),rstack(irstack),istack(iistack))
          do iclv=1,2
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('eclvp')    ! customise Eirene contour level by proportion
          stack_req=(/2,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          ix=ipop(istack,iistack)
          iy=ipop(istack,iistack)
          call eclvp(ix,iy)
          write(stmp,'(a,2x,a1,2i4,a2)')
     .     trim(eirene_label(iewrk)),'{',ix,iy,' }'
          eirene_label(iewrk)=stmp
        case('eabs')     ! abs() on Eirene stack
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,nedim(2,iewrk)-1
            do i=-1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=abs(eirene_stack(i,is,iewrk))
            end do
          end do
          eirene_label(iewrk)=trim(eirene_label(iewrk))//' (ABS)'
        case('esumz')    ! sum over ns for Eirene stack data
          stack_req=(/2,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          izstrt=istack(iistack-1)
          izend=istack(iistack)
          if(izstrt .eq. 0 .and. izend .eq. 0) then
            izstrt=0
            izend=nedim(2,iewrk)-1
          end if
          if(izend .eq. -1) izend=nedim(2,iewrk)-1
          if(izstrt .lt. 0 .or. izend .gt. nedim(2,iewrk)-1 .or.
     1       izstrt .gt. izend) then
            write(*,*) 'Parameters supplied to esumz were invalid - ',
     1                 'valid parameters are:'
            write(*,*) ' 0 <= ',izstrt,' <= ',izend,' <= ',
     1           nedim(2,iewrk)-1
            return
          end if
          do is=izstrt+1,izend
            do i=1,nedim(1,iewrk)
              eirene_stack(i,izstrt,iewrk)=
     1          eirene_stack(i,izstrt,iewrk)
     2         +eirene_stack(i,is,iewrk)
            end do
          end do
          do i=1,nedim(1,iewrk)
            eirene_stack(i,0,iewrk)=
     1       eirene_stack(i,izstrt,iewrk)
          end do
          nedim(2,iewrk)=1
          if(izstrt .ne. izend) then
            eirene_label(iewrk)=
     1       trim(eirene_label(iewrk))
     2       //' SUM'//
     3       trim(textewk(izstrt,iewrk))
     4       //' to'//
     5       trim(textewk(izend,iewrk))
          else
            wall_label(iewrk)=
     1       trim(wall_label(iewrk)) // trim(textewk(izstrt,iewrk))
          end if
          textewk(0,iewrk)=textbl(0)
          call pop_stack(stack_req)
        case('em+')      ! Eirene stack addition
          stack_req=(/0,0,0,0,2,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          if(nedim(2,iewrk-1) .ne. nedim(2,iewrk) .and.
     1       nedim(2,iewrk-1) .ne. 1 .and. nedim(2,iewrk) .ne. 1) then
            write(*,*) 'Mismatch in nedim(2)''s in operands'
            goto 4900
          end if
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk+1)=
     1         eirene_stack(i,min(is,nedim(2,iewrk-1)-1),iewrk-1)+
     2         eirene_stack(i,min(is,nedim(2,iewrk)-1),iewrk)
            end do
          end do
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk-1)=
     1         eirene_stack(i,is,iewrk+1)
            end do
            if(nedim(2,iewrk-1) .gt. nedim(2,iewrk)) then
              textewk(is,iewrk-1)=textewk(is,iewrk-1)
            else
              textewk(is,iewrk-1)=textewk(is,iewrk)
            end if
          end do
          eirene_label(iewrk-1)='Special'
          nedim(2,iewrk-1)=max(nedim(2,iewrk-1),nedim(2,iewrk))
          iewrk=iewrk-1
        case('em-')      ! Eirene stack subtraction
          stack_req=(/0,0,0,0,2,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          if(nedim(2,iewrk-1) .ne. nedim(2,iewrk) .and.
     1       nedim(2,iewrk-1) .ne. 1 .and. nedim(2,iewrk) .ne. 1) then
            write(*,*) 'Mismatch in nedim(2)''s in operands'
            goto 4900
          end if
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk+1)=
     1         eirene_stack(i,min(is,nedim(2,iewrk-1)-1),iewrk-1)-
     2         eirene_stack(i,min(is,nedim(2,iewrk)-1),iewrk)
            end do
          end do
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk-1)=eirene_stack(i,is,iewrk+1)
            end do
            if(nedim(2,iewrk-1) .gt. nedim(2,iewrk)) then
              textewk(is,iewrk-1)=textewk(is,iewrk-1)
            else
              textewk(is,iewrk-1)=textewk(is,iewrk)
            end if
          end do
          eirene_label(iewrk-1)='Special'
          nedim(2,iewrk-1)=max(nedim(2,iewrk-1),nedim(2,iewrk))
          iewrk=iewrk-1
        case('em*')      ! Eirene matrix multiplication
          stack_req=(/0,0,0,0,2,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          if(nedim(2,iewrk-1) .ne. nedim(2,iewrk) .and.
     1       nedim(2,iewrk-1) .ne. 1 .and. nedim(2,iewrk) .ne. 1) then
            write(*,*) 'Mismatch in nedim(2)''s in operands'
            return
          end if
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk+1)=
     1         eirene_stack(i,min(is,nedim(2,iewrk-1)-1),iewrk-1)*
     2         eirene_stack(i,min(is,nedim(2,iewrk)-1),iewrk)
            end do
          end do
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk-1)=eirene_stack(i,is,iewrk+1)
            end do
            if(nedim(2,iewrk-1) .gt. nedim(2,iewrk)) then
              textewk(is,iewrk-1)=textewk(is,iewrk-1)
            else
              textewk(is,iewrk-1)=textewk(is,iewrk)
            end if
          end do
          eirene_label(iewrk-1)='Special'
          nedim(2,iewrk-1)=max(nedim(2,iewrk-1),nedim(2,iewrk))
          iewrk=iewrk-1
        case('em/')      ! Eirene matrix division
          stack_req=(/0,0,0,0,2,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          if(nedim(2,iewrk-1) .ne. nedim(2,iewrk) .and.
     1       nedim(2,iewrk-1) .ne. 1 .and. nedim(2,iewrk) .ne. 1) then
            write(*,*) 'Mismatch in nedim(2)''s in operands'
            return
          end if
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i=1,nedim(1,iewrk)
              if (eirene_stack(i,
     .             min(is,nedim(2,iewrk)-1),iewrk).ne.0.0_R8) then
                eirene_stack(i,is,iewrk+1)=
     1           eirene_stack(i,min(is,nedim(2,iewrk-1)-1),iewrk-1)/
     2           eirene_stack(i,min(is,nedim(2,iewrk)-1),iewrk)
              else
                eirene_stack(i,is,iewrk+1)=0.0_R8
              end if
            end do
          end do
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk-1)=eirene_stack(i,is,iewrk+1)
            end do
            if(nedim(2,iewrk-1) .gt. nedim(2,iewrk)) then
              textewk(is,iewrk-1)=textewk(is,iewrk-1)
            else
              textewk(is,iewrk-1)=textewk(is,iewrk)
            end if
          end do
          eirene_label(iewrk-1)='Special'
          nedim(2,iewrk-1)=max(nedim(2,iewrk-1),nedim(2,iewrk))
          iewrk=iewrk-1
        case('rem+')     ! Eirene stack + real
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          do is=0,nedim(2,iewrk)-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=
     1         eirene_stack(i,is,iewrk)+rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('rem-')     ! Eirene stack - real
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          do is=0,nedim(2,iewrk)-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=
     1         eirene_stack(i,is,iewrk)-rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('rem*')     ! Eirene stack * real
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          do is=0,nedim(2,iewrk)-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=
     1         eirene_stack(i,is,iewrk)*rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('rem/')     ! Eirene stack / real
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          if(rstack(irstack).eq.0.0_R8) then
            write (*,*) 'Divide by zero, You lose'
            goto 4900
          endif
          do is=0,nedim(2,iewrk)-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=
     1         eirene_stack(i,is,iewrk)/rstack(irstack)
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('rem**')    ! Eirene stack ** real
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 5900
          do is=0,nedim(2,iewrk)-1
            do i=1,nedim(1,iewrk)
              if (rstack(irstack).ne.0.0_R8) then
                eirene_stack(i,is,iewrk)=
     1           eirene_stack(i,is,iewrk)**rstack(irstack)
              else
                eirene_stack(i,is,iewrk)=1.0_R8
              end if
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('edrop')    ! drop top of Eirene stack
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          iewrk=iewrk-1
        case('edup')     ! duplicate top of Eirene stack
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          if(iewrk .ge. wklng) call dropeirenestackbottom()
          iewrk=iewrk+1
          nedim(1,iewrk)=nedim(1,iewrk-1)
          nedim(2,iewrk)=nedim(2,iewrk-1)
          do is=0,nedim(2,iewrk-1)-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=eirene_stack(i,is,iewrk-1)
            end do
            textewk(is,iewrk)=textewk(is,iewrk-1)
          end do
          eirene_label(iewrk)=eirene_label(iewrk-1)
        case('eswap')    ! swap top two elements of Eirene stack
          stack_req=(/0,0,0,0,2,0,0/)
          if(bad_args(stack_req)) return   ! 4900
          call iswap(nedim(1,iewrk),nedim(1,iewrk-1))
          call iswap(nedim(2,iewrk),nedim(2,iewrk-1))
          do is=0,nedim(2,iewrk-1)-1
            do i=1,nedim(1,iewrk)
              call rswap(eirene_stack(i,is,iewrk),
     1                   eirene_stack(i,is,iewrk-1))
            end do
            call strswap(textewk(is,iewrk),textewk(is,iewrk-1))
          end do
          call strswap(eirene_label(iewrk),eirene_label(iewrk-1))
        case('eover')    ! e2 e1 --> e2 e1 e2
          stack_req=(/0,0,0,0,2,0,0/)
          if(bad_args(stack_req)) return ! 4900
          if(iewrk .ge. wklng) call dropeirenestackbottom()
          iewrk=iewrk+1
          nedim(1,iewrk)=nedim(1,iewrk-2)
          nedim(2,iewrk)=nedim(2,iewrk-2)
          do is=0,nedim(2,iewrk-2)-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=eirene_stack(i,is,iewrk-2)
            end do
            textewk(is,iewrk)=textewk(is,iewrk-2)
          end do
          eirene_label(iewrk)=eirene_label(iewrk-2)
        case('elist')    ! list contents of Eirene stack
          if(iewrk.eq.0) then
            write(*,*) 'No data on Eirene stack'
            return
          end if
          do i = iewrk, 1, -1
            write(*,*) trim(eirene_label(i)),
     &       (trim(textewk(is,i)),is=0,nedim(2,i)-1)
          end do
        case('esel')     ! select Eirene element
          stack_req=(/1,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          i1=istack(iistack)
          if(1.gt.i1.or.i1.gt.nedim(1,iewrk)) then
            write(*,*) 'Faulty parameter i1 = ',i1,' for esel!'
            return
          endif
          do is=0,nedim(2,iewrk)-1
            do i=1,nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=eirene_stack(i1,is,iewrk)
            end do
          end do
          call pop_stack(stack_req)
        case('ezsel')    ! select species in Eirene stack
          stack_req=(/1,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          is1=istack(iistack)
          if(0.gt.is1.or.is1.ge.nedim(2,iewrk)) then
            write(*,*) 'Faulty parameter is = ',is1,' for ezsel!'
            return
          endif
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=eirene_stack(i,is1,iewrk)
              earrowx(i,is)=earrowx(i,is1)
              earrowy(i,is)=earrowy(i,is1)
            end do
            textewk(is,iewrk)=textewk(is1,iewrk)
          end do
          nedim(2,iewrk)=1
          call pop_stack(stack_req)
        case('emmin')    ! min(e1,e2)
          stack_req=(/0,0,0,0,2,0,0/)
          if(bad_args(stack_req)) return ! 4900
          if(nedim(2,iewrk-1) .ne. nedim(2,iewrk) .and.
     1       nedim(2,iewrk-1) .ne. 1 .and. nedim(2,iewrk) .ne. 1) then
            write(*,*) 'Mismatch in nedim(2)''s in operands'
            goto 4900
          end if
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i = 1, nedim(1,iewrk)
              eirene_stack(i,is,iewrk+1)=
     1         min(eirene_stack(i,min(is,nedim(2,iewrk-1)-1),iewrk-1),
     2             eirene_stack(i,min(is,nedim(2,iewrk)-1),iewrk))
            end do
          end do
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i = 1, nedim(1,iewrk)
              eirene_stack(i,is,iewrk-1)=eirene_stack(i,is,iewrk+1)
            end do
            if(nedim(2,iewrk-1) .gt. nedim(2,iewrk)) then
              textewk(is,iewrk-1)=textewk(is,iewrk-1)
            else
              textewk(is,iewrk-1)=textewk(is,iewrk)
            end if
          end do
          eirene_label(iewrk-1)='Special'
          nedim(2,iewrk-1)=max(nedim(2,iewrk-1),nedim(2,iewrk))
          iewrk=iewrk-1
        case('emmax')    ! max(e1,e2)
          stack_req=(/0,0,0,0,2,0,0/)
          if(bad_args(stack_req)) return ! 4900
          if(nedim(2,iewrk-1) .ne. nedim(2,iewrk) .and.
     1       nedim(2,iewrk-1) .ne. 1 .and. nedim(2,iewrk) .ne. 1) then
            write(*,*) 'Mismatch in nedim(2)''s in operands'
            goto 4900
          end if
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i = 1, nedim(1,iewrk)
              eirene_stack(i,is,iewrk+1)=
     1         max(eirene_stack(i,min(is,nedim(2,iewrk-1)-1),iewrk-1),
     2             eirene_stack(i,min(is,nedim(2,iewrk)-1),iewrk))
            end do
          end do
          do is=0,max(nedim(2,iewrk-1),nedim(2,iewrk))-1
            do i = 1, nedim(1,iewrk)
              eirene_stack(i,is,iewrk-1)=eirene_stack(i,is,iewrk+1)
            end do
            if(nedim(2,iewrk-1) .gt. nedim(2,iewrk)) then
              textewk(is,iewrk-1)=textewk(is,iewrk-1)
            else
              textewk(is,iewrk-1)=textewk(is,iewrk)
            end if
          end do
          eirene_label(iewrk-1)='Special'
          nedim(2,iewrk-1)=max(nedim(2,iewrk-1),nedim(2,iewrk))
          iewrk=iewrk-1
        case('remmi')    ! min(r1,e1)
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 4900
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=
     1         min(eirene_stack(i,is,iewrk),rstack(irstack))
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('remma')    ! max(r1,e1)
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 4900
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=
     1         max(eirene_stack(i,is,iewrk),rstack(irstack))
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('emlt')     ! e<r
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 5900
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              if(eirene_stack(i,is,iewrk).lt.
     1                 rstack(irstack)) then
                eirene_stack(i,is,iewrk)=1.0_R8
              else
                eirene_stack(i,is,iewrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('emle')     ! e<=r
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 4900
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              if(eirene_stack(i,is,iewrk).le.
     1                 rstack(irstack)) then
                eirene_stack(i,is,iewrk)=1.0_R8
              else
                eirene_stack(i,is,iewrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('emeq')     ! e=r
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 4900
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              if(eirene_stack(i,is,iewrk).eq.
     1                 rstack(irstack)) then
                eirene_stack(i,is,iewrk)=1.0_R8
              else
                eirene_stack(i,is,iewrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('emne')     ! e!=r
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 4900
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              if(eirene_stack(i,is,iewrk).ne.
     1                 rstack(irstack)) then
                eirene_stack(i,is,iewrk)=1.0_R8
              else
                eirene_stack(i,is,iewrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('emgt')     ! e>r
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 5900
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              if(eirene_stack(i,is,iewrk).gt.
     1                 rstack(irstack)) then
                eirene_stack(i,is,iewrk)=1.0_R8
              else
                eirene_stack(i,is,iewrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('emge')     ! e=>r
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 4900
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              if(eirene_stack(i,is,iewrk).ge.
     1                 rstack(irstack)) then
                eirene_stack(i,is,iewrk)=1.0_R8
              else
                eirene_stack(i,is,iewrk)=0.0_R8
              endif
            end do
          end do
          dummy=rpop(rstack,irstack)
          eirene_label(iewrk)='Special'
        case('elog')     ! log10(Eirene)
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              if(eirene_stack(i,is,iewrk).gt.0.0_R8) then
                eirene_stack(i,is,iewrk)=log10(eirene_stack(i,is,iewrk))
              else
                write(*,*) 'Negative argument for ',
     1            i,is,iewrk,' : Nothing done.'
              end if
            end do
          end do
          eirene_label(iewrk)='Special'
        case('eln')      ! ln(Eirene)
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              if(eirene_stack(i,is,iewrk).gt.0.0_R8) then
                eirene_stack(i,is,iewrk)=log(eirene_stack(i,is,iewrk))
              else
                write(*,*) 'Negative argument for ',
     1           i,is,iewrk,'. Nothing done.'
              end if
            end do
          end do
          eirene_label(iewrk)='Special'
        case('esqrt')    ! sqrt(Eirene)
           stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              if(eirene_stack(i,is,iewrk).ge.0.0_R8) then
                eirene_stack(i,is,iewrk)=sqrt(eirene_stack(i,is,iewrk))
              else
                write(*,*) 'Negative argument for ',
     1            i,is,iewrk,' : Nothing done.'
              end if
            end do
          end do
          eirene_label(iewrk)='Special'
        case('eexp')     ! exp(Eirene)
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=exp(eirene_stack(i,is,iewrk))
            end do
          end do
          eirene_label(iewrk)='Special'
        case('ezrng')    ! select species range in Eirene stack
          stack_req=(/2,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          if(istack(iistack-1).gt.istack(iistack)) then
            write(*,*) 'Parameters for ezrng must be given in ',
     &       'ascending order!'
            return
          end if
          if(0.gt.istack(iistack-1).or.
     &       istack(iistack-1).ge.nedim(2,iewrk)) then
            write(*,*) 'Faulty first parameter is = ',istack(iistack-1),
     &       ' for ezrng!'
            return
          endif
          if(0.gt.istack(iistack).or.
     &       istack(iistack).ge.nedim(2,iewrk)) then
            write(*,*) 'Faulty second parameter is = ',istack(iistack),
     &       ' for ezrng!'
            return
          endif
          do is=istack(iistack-1),istack(iistack)
            do i = 1, nedim(1,iewrk)
              eirene_stack(i,is-istack(iistack-1),iewrk)=
     1         eirene_stack(i,is,iewrk)
              earrowx(i,is-istack(iistack-1))=earrowx(i,is)
              earrowy(i,is-istack(iistack-1))=earrowy(i,is)
            end do
            textewk(is-istack(iistack-1),iewrk)=textewk(is,iewrk)
          end do
          nedim(2,iewrk)=istack(iistack)-istack(iistack-1)+1
          call pop_stack(stack_req)
        case('echs')     ! change sign in Eirene stack
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,nedim(2,iewrk)-1
            do i = 1, nedim(1,iewrk)
              eirene_stack(i,is,iewrk)=-eirene_stack(i,is,iewrk)
            end do
          end do
          eirene_label(iewrk)=trim(eirene_label(iewrk))//' (CHS)'
        case('etos')     ! Eirene array -> m
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          if(iwrk .eq. wklng) call dropstackbottom()
          iwrk=iwrk+1
          ndim(1,iwrk)=nx
          ndim(2,iwrk)=ny
          ndim(3,iwrk)=nedim(2,iewrk)
          work(-1:nx,-1:ny,0:ndim(3,iwrk)-1,iwrk)=0.0_R8
          nntria(0:nx-1,0:ny-1)=0
          do i=1,ntrii
            if (b2_cell(ixtri(i),iytri(i))) then
              nntria(ix_e2b(ixtri(i)),iytri(i)-1)=
     .         nntria(ix_e2b(ixtri(i)),iytri(i)-1)+1
              work(ix_e2b(ixtri(i)),iytri(i)-1,0:ndim(3,iwrk)-1,iwrk)=
     .         work(ix_e2b(ixtri(i)),iytri(i)-1,0:ndim(3,iwrk)-1,iwrk)+
     .         eirene_stack(i,0:ndim(3,iwrk)-1,iewrk)
            end if
          enddo
          do iy = 0, ny-1
            do ix = 0, nx-1
              if (nntria(ix,iy).gt.0) then
                work(ix,iy,0:ndim(3,iwrk)-1,iwrk)=
     1           work(ix,iy,0:ndim(3,iwrk)-1,iwrk)/nntria(ix,iy)
              end if
            end do
          end do
          call extend(work(-1:nx,-1:ny,0:ndim(3,iwrk)-1,iwrk),
     1     nx,ny,ndim(3,iwrk),nx,ny,ndim(3,iwrk))
          values(9)=minval(work(-1:ndim(1,iwrk),-1:ndim(2,iwrk),
     1         0:ndim(3,iwrk)-1,iwrk))
          values(10)=maxval(work(-1:ndim(1,iwrk),-1:ndim(2,iwrk),
     1         0:ndim(3,iwrk)-1,iwrk))
          textwk(0:ndim(3,iwrk)-1,iwrk)=textewk(0:ndim(3,iwrk)-1,iewrk)
          graphlabel(iwrk)=eirene_label(iewrk)
          goto 4049
        case('epos')     ! Eirene stack with 1 if element > 0, else 0
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          do is=0,nedim(2,iewrk)-1
            do i=1,nedim(1,iewrk)
              if(eirene_stack(i,is,iewrk) .gt. 0.0_R8) then
                eirene_stack(i,is,iewrk)=1.0_R8
              else
                eirene_stack(i,is,iewrk)=0.0_R8
              endif
            end do
          end do
          eirene_label(iewrk)='Special'
        case('evac')     ! fill vacuum region of Eirene stack with a value
          stack_req=(/0,1,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 4900
          do is=0,nedim(2,iewrk)-1
            do i=1,nedim(1,iewrk)
              if (.not.b2_cell(ixtri(i),iytri(i))) then
                eirene_stack(i,is,iewrk) = rstack(irstack)
              end if
            enddo
          enddo
          dummy=rpop(rstack,irstack)
        case('earox')    ! set arrow field x-component from Eirene data
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 4900
          do is=0,max(natmi,nmoli,nioni)-1
            do i=1,nedim(1,iewrk)
              earrowx(i,is)=
     1         eirene_stack(i,min(is,nedim(2,iewrk)-1),iewrk)
            enddo
          enddo
          goto 4049
        case('earoy')    ! set arrow field y-component from Eirene data
          stack_req=(/0,0,0,0,1,0,0/)
          if(bad_args(stack_req)) return ! 4900
          do is=0,max(natmi,nmoli,nioni)-1
            do i=1,nedim(1,iewrk)
              earrowy(i,is)=
     1         eirene_stack(i,min(is,nedim(2,iewrk)-1),iewrk)
            enddo
          enddo
          goto 4049
        case('etal')     ! read Eirene tally output file onto the Eirene stack
          stack_req=(/0,0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          arg=spop(sstack,isstack)
          call find_file(arg,ex)
          if (.not.ex) then
            write(*,*) 'Eirene tally file ', trim(arg), ' not found!'
            return
          end if
          nz=0
          open(99,file=arg,status='old')
          iewrk=iewrk+1
          if(iewrk .eq. wklng) call dropeirenestackbottom()
          nedim(1,iewrk)=ntrii
 4001     continue
          read(99,1028,end=4003)
          read(99,1028)
          read(99,'(a)') stmp  ! Tally title
          read(99,'(a)') stmp1 ! Species name
          read(99,'(a)') stmp2 ! Units
          read(99,1028)
          read(99,1028)
          read(99,1028)
          call strip_spaces(stmp2)
          description=trim(stmp)//' ('//trim(stmp2)//')'
          textewk(nz,iewrk)=trim(stmp1)
          do i=1,ntrii
            read(99,*) j,dummy,eirene_stack(j,nz,iewrk)
          end do
          read(99,1028)
          read(99,1028) ! average value
          read(99,1028)
          read(99,1028) ! additional cells
 4002     continue
          read(99,'(a)',end=4003) stmp
          if (index(stmp,'=').eq.0) then
            goto 4002
          else
            nz=nz+1
          end if
          goto 4001
 4003     continue
          nedim(2,iewrk)=nz
          values(38)=minval(eirene_stack(1:ntrii,0:nz-1,iewrk))
          values(39)=minval(eirene_stack(1:ntrii,0:nz-1,iewrk))
          call customise_label(description,eirene_label(iewrk))
          call set_wallbarlabel(eirene_label(iewrk))
          if (b2plot_debug.ge.1)
     .     write(*,'(4a)') 'Copied ',trim(description),
     .                     ' to Eirene stack from file ',trim(arg)
        case('echor')    ! chord file for Eirene stack data
          stack_req=(/0,0,1,0,1,0,0/)
          if(bad_args(stack_req)) return
          call process_eirene_chords(spop(sstack,isstack))
        case('echbx')    ! boxcar average of chords for Eirene stack data
          stack_req=(/0,0,1,0,1,0,0/)
          if(bad_args(stack_req)) return
          call process_eirene_chords_boxcar(spop(sstack,isstack))
        case('echvl')    ! Eirene stack values at positions
          stack_req=(/0,0,1,0,1,0,0/)
          if(bad_args(stack_req)) return
          call process_eirene_chord_slice(spop(sstack,isstack))
        case('echmp')    ! plots a contour plot of the chord contribution to the Eirene element
          stack_req=(/0,0,1,0,1,0,0/)
          if(bad_args(stack_req)) return
          call map_eirene_chords(spop(sstack,isstack))
        case('ecam')     ! camera plot of Eirene stack data
          stack_req=(/0,0,1,0,1,0,0/)
          if(bad_args(stack_req)) return
          call process_eirene_camera(spop(sstack,isstack))
        case('eoutp')    ! prints the value of a specific Eirene element
          stack_req=(/1,0,1,0,1,0,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          it=ipop(istack,iistack)
          if(0.gt.it.or.it.gt.nedim(1,iewrk)) then
            write(*,*) 'Faulty parameter it = ',it,' for eoutp!'
            return
          endif
          do is=0,nedim(2,iewrk)-1
            if(nedim(2,iewrk) .eq. 1) then
              write(*,'(a,1x,1p,g15.8)')
     1                 trim(stmp)//trim(textewk(is,iewrk)),
     1                 eirene_stack(it,is,iewrk)
              if(writing_to_file) then
                write(nwtf,'(a,1x,1p,g15.8)')
     1                 trim(stmp)//trim(textewk(is,iewrk)),
     1                 eirene_stack(it,is,iewrk)
              endif
            else
              write(*,'(a,i3,1p,g15.8)')
     1               trim(stmp)//trim(textewk(is,iewrk)),is,
     1               eirene_stack(it,is,iewrk)
              if(writing_to_file) then
                write(nwtf,'(a,i3,1p,g15.8)')
     1               trim(stmp)//trim(textewk(is,iewrk)),is,
     1               eirene_stack(it,is,iewrk)
              endif
            endif
          enddo
        case('eprt')     ! output the top of Eirene stack to a file
          stack_req=(/0,0,1,0,1,0,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          open(99,file=stmp)
          do is=0,nedim(2,iewrk)-1
            write(99,'(''#'',1x,a,1x,a)')
     1       trim(eirene_label(iewrk)),trim(textewk(is,iewrk))
            do it=1,ntrii
              write(99,*) it,is,eirene_stack(it,is,iewrk)
            enddo
          enddo
          close(99)
        case('eavg')
          stack_req=(/0,8,0,0,1,0,0/)
          if(bad_args(stack_req)) return
          volsum=0.0_R8
          allocate(eavg(0:nedim(2,iewrk)-1))
          eavg=0.0_R8
          y4 = rpop(rstack,irstack)
          x4 = rpop(rstack,irstack)
          y3 = rpop(rstack,irstack)
          x3 = rpop(rstack,irstack)
          y2 = rpop(rstack,irstack)
          x2 = rpop(rstack,irstack)
          y1 = rpop(rstack,irstack)
          x1 = rpop(rstack,irstack)
          do i=1,ntrii
            xc=(xtrian(necke(1,i))+xtrian(necke(2,i))+
     .          xtrian(necke(3,i)))/3.0_R8
            yc=(ytrian(necke(1,i))+ytrian(necke(2,i))+
     .          ytrian(necke(3,i)))/3.0_R8
            if (inside_rectangle(xc,yc,x1,y1,x2,y2,x3,y3,x4,y4)) then
              volsum=volsum+triangle_vol(i)
              do is=0,nedim(2,iewrk)-1
                eavg(is)=eavg(is)+
     .           eirene_stack(i,is,iewrk)*triangle_vol(i)
              end do
            end if
          end do
          if (volsum.gt.0.0_R8) then
            eavg(:)=eavg(:)/volsum
            write(6,'(3a)')
     .       'Volume-averaged ',trim(eirene_label(iewrk)),' in region '
            write(6,'(1p,a,8(e10.3,a))')
     .       ' [(',x1,',',y1,');(',x2,',',y2,');(',
     .             x3,',',y3,');(',x4,',',y4,')] : '
            do is = 0, nedim(2,iewrk)-1
              write(6,'(a,1p,1x,e10.3)') textewk(is,iewrk),eavg(is)
            end do
          else
            write(6,*) 'No triangle centroids found in chosen region!'
          end if
          deallocate(eavg)
        case default EFVALUE
          cflag=.false.
      END SELECT EFVALUE
c
      goto 4900
 1028 format()

 4049 continue
      stack_req=(/0,0,0,0,1,0,0/)
      if(bad_args(stack_req)) return
      iewrk=iewrk-1

 4900 continue
      if(cflag) then
        if(b2plot_debug.ge.4)
     1   write(*,*) 'Applied operator ',trim(cval),
     2   ' to Eirene stack array'
        if(iewrk.gt.0) then
          values(38)=minval(eirene_stack(1:nedim(1,iewrk),
     1         0:nedim(2,iewrk)-1,iewrk))
          values(39)=maxval(eirene_stack(1:nedim(1,iewrk),
     1         0:nedim(2,iewrk)-1,iewrk))
        endif
      endif
#endif
      return

      end subroutine select_eirene_command

      subroutine select_wall_stack
      integer i, ibnd, itrack

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
          description='Wall particle flux (s^{-1})'
          nz=nsw
          call wall_copy(wall_flux,1.0_R8,nz,textwl)
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
          stack_req=(/1,0,0,0,0,0,0/)
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
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          data(1:istack(iistack)-1) = 1.0_R8
          data(istack(iistack):nbnd) = 0.0_R8
          write(*,*) '"iwlt" option with ',ipop(istack,iistack)
          call wall_copy(data,1.0_R8,nz,textbl)
        case ('iwle')
          description='iwle'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          data(1:istack(iistack)) = 1.0_R8
          data(istack(iistack)+1:nbnd) = 0.0_R8
          write(*,*) '"iwle" option with ',ipop(istack,iistack)
          call wall_copy(data,1.0_R8,nz,textbl)
        case ('iweq')
          description='iweq'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          data(1:nbnd) = 0.0_R8
          data(istack(iistack)) = 1.0_R8
          write(*,*) '"iweq" option with ',ipop(istack,iistack)
          call wall_copy(data,1.0_R8,nz,textbl)
        case ('iwne')
          description='iwne'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          data(1:nbnd) = 1.0_R8
          data(istack(iistack)) = 0.0_R8
          write(*,*) '"iwne" option with ',ipop(istack,iistack)
          call wall_copy(data,1.0_R8,nz,textbl)
        case ('iwgt')
          description='iwgt'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          data(1:istack(iistack)) = 0.0_R8
          data(istack(iistack)+1:nbnd) = 1.0_R8
          write(*,*) '"iwgt" option with ',ipop(istack,iistack)
          call wall_copy(data,1.0_R8,nz,textbl)
        case ('iwge')
          description='iwge'
          nz=1
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          data(1:istack(iistack)-1) = 0.0_R8
          data(istack(iistack):nbnd) = 1.0_R8
          write(*,*) '"iwge" option with ',ipop(istack,iistack)
          call wall_copy(data,1.0_R8,nz,textbl)
        case ('iswlt')
          description='iswlt'
          nz=nsw
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wall_tmp(1:nbnd,0:istack(iistack)-1) = 1.0_R8
          wall_tmp(1:nbnd,istack(iistack):nz-1) = 0.0_R8
          write(*,*) '"iswlt" option with ',ipop(istack,iistack)
          call wall_copy(wall_tmp,1.0_R8,nz,textwl)
        case ('iswle')
          description='iswle'
          nz=nsw
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wall_tmp(1:nbnd,0:istack(iistack)) = 1.0_R8
          wall_tmp(1:nbnd,istack(iistack)+1:nz-1) = 0.0_R8
          write(*,*) '"iswle" option with ',ipop(istack,iistack)
          call wall_copy(wall_tmp,1.0_R8,nz,textwl)
        case ('isweq')
          description='isweq'
          nz=nsw
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wall_tmp(1:nbnd,0:nz-1) = 0.0_R8
          wall_tmp(1:nbnd,istack(iistack)) = 1.0_R8
          write(*,*) '"isweq" option with ',ipop(istack,iistack)
          call wall_copy(wall_tmp,1.0_R8,nz,textwl)
        case ('iswne')
          description='iswne'
          nz=nsw
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wall_tmp(1:nbnd,0:nz-1) = 1.0_R8
          wall_tmp(1:nbnd,istack(iistack)) = 0.0_R8
          write(*,*) '"iswne" option with ',ipop(istack,iistack)
          call wall_copy(wall_tmp,1.0_R8,nz,textwl)
        case ('iswgt')
          description='iswgt'
          nz=nsw
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wall_tmp(1:nbnd,0:istack(iistack)) = 0.0_R8
          wall_tmp(1:nbnd,istack(iistack)+1:nz-1) = 1.0_R8
          write(*,*) '"iswgt" option with ',ipop(istack,iistack)
          call wall_copy(wall_tmp,1.0_R8,nz,textwl)
        case ('iswge')
          description='iswge'
          nz=nsw
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wall_tmp(1:nbnd,0:istack(iistack)-1) = 0.0_R8
          wall_tmp(1:nbnd,istack(iistack):nz-1) = 1.0_R8
          write(*,*) '"iswge" option with ',ipop(istack,iistack)
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
          description='Physical sputtering yield'
          if (sput_frac_flag.eq.0) then
            if (sput_dst.lt.0 .or. sput_dst.ge.ns) return
            description=trim(description)//' of'//
     .       trim(textpl(sput_dst))//' by'
          else
            description=trim(description)//' of wall materials by'
          endif
          nz=nsw
          call wall_copy(physical_sputtering,1.0_R8,nz,textwl)
        case('sput2')
          description='Physical sputtering energy yield'
          if (sput_frac_flag.eq.0) then
            if (sput_dst.lt.0 .or. sput_dst.ge.ns) return
            description=trim(description)//' of' //
     .       trim(textpl(sput_dst)) //' by'
          else
            description=trim(description)//' of wall materials by'
          endif
          nz=nsw
          call wall_copy(physical_sputtering_energy,1.0_R8,nz,textwl)
        case('sputchem')
          description='Chemical sputtering yield'
          nz=nsw
          call wall_copy(chemical_sputtering,1.0_R8,nz,textwl)
        case('thermev')
          description='Thermal evaporation fluxes (/m^{2}/s)'
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
          call wall_copy(erosion,1.0_R8,nz,texttargsp)
        case('deposit')
          description='Gross deposition'
          nz=ntargsp
          call wall_copy(deposition,1.0_R8,nz,texttargsp)
        case('ero_ml')
          description='Monolayer erosion'
          nz=ntargsp
          call wall_copy(monolayer_erosion,1.0_R8,nz,texttargsp)
        case('depo_ml')
          description='Monolayer deposition'
          nz=ntargsp
          call wall_copy(monolayer_deposition,1.0_R8,nz,texttargsp)
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
        case('kappa_a5')
          description='Bulk Kappa A5 polynomial factor (W/m/K^{5})'
          nz=1
          call wall_copy(matbulk_Kappa_A5,1.0_R8,nz,textbl)
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
        case('kappapa5')
          description='Bulk Kappa_perp A5 polynomial factor (W/m/K^{5})'
          nz=1
          call wall_copy(matbulk_Kappa_perp_A5,1.0_R8,nz,textbl)
        case('tgdns')
          description='Bulk mass density (kg/m^{3})'
          nz=1
          call wall_copy(matbulk_target_density,1.0_R8,nz,textbl)
        case('cp_b1')
          description='Bulk heat capacity B1 factor (J/kg/K)'
          nz=1
          call wall_copy(matbulk_Cp_B1,1.0_R8,nz,textbl)
        case('cp_b2')
          description='Bulk heat capacity B2 factor (J/kg/K^{2})'
          nz=1
          call wall_copy(matbulk_Cp_B2,1.0_R8,nz,textbl)
        case('cp_b3')
          description='Bulk heat capacity B3 factor'
          if (maxval(matbulk_format(1:nbnd)).eq.0) then
            description=trim(description)//' (J.K/kg)'
          else if (minval(matbulk_format(1:nbnd)).eq.1) then
            description=trim(description)//' (J/kg/K^{3})'
          else
            write(*,*) 'Warning: inconsistent units for cp_b3 !'
          endif
          nz=1
          call wall_copy(matbulk_Cp_B3,1.0_R8,nz,textbl)
        case('cp_b4')
          description='Bulk heat capacity B4 factor (J/kg/K^{4})'
          nz=1
          call wall_copy(matbulk_Cp_B4,1.0_R8,nz,textbl)
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
        case('ckappaa5')
          description='Coating Kappa A5 polynomial factor (W/m/K^{5})'
          nz=1
          call wall_copy(coating_Kappa_A5,1.0_R8,nz,textbl)
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
        case('ckapppa5')
          description=
     .     'Coating Kappa_perp A5 polynomial factor (W/m/K^{5})'
          nz=1
          call wall_copy(coating_Kappa_perp_A5,1.0_R8,nz,textbl)
        case('ctgdns')
          description='Coating mass density (kg/m^{3})'
          nz=1
          call wall_copy(coating_target_density,1.0_R8,nz,textbl)
        case('c_cp_b1')
          description='Coating heat capacity B1 factor (J/kg/K)'
          nz=1
          call wall_copy(coating_Cp_B1,1.0_R8,nz,textbl)
        case('c_cp_b2')
          description='Coating heat capacity B2 factor (J/kg/K^{2})'
          nz=1
          call wall_copy(coating_Cp_B2,1.0_R8,nz,textbl)
        case('c_cp_b3')
          description='Coating heat capacity B3 factor'
          if (maxval(coating_format(1:nbnd)).eq.0) then
            description=trim(description)//' (J.K/kg)'
          else if (minval(coating_format(1:nbnd)).eq.1) then
            description=trim(description)//' (J/kg/K^{3})'
          else
            write(*,*) 'Warning: inconsistent units for c_cp_b3 !'
          endif
          nz=1
          call wall_copy(coating_Cp_B3,1.0_R8,nz,textbl)
        case('c_cp_b4')
          description='Coating heat capacity B4 factor (J/kg/K^{4})'
          nz=1
          call wall_copy(coating_Cp_B4,1.0_R8,nz,textbl)
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
            data(ibnd)=qrad(ibnd,target_temp(ibnd,1),.true.)
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
      integer i, is, ix, iclv, iseg

      cflag=.true.
      WFVALUE: select case (cval)
        case('wall')     ! adds wall stack data to contour plot
          plotting_wall=.true.
          call set_wall_on(values(14),values(15),wlinlog)
        case('wall_ank') ! add Kukushkin wall data to contour plot
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          plotting_wall=istack(iistack) .ne. 0
          call wallon(.true.,istack(iistack),values(14),values(15),
     .     wlinlog)
          call pop_stack(stack_req)
        case('logw')     ! log wall elements
          if (mod(wlinlog,2).eq.1) wlinlog=wlinlog+1
          call set_wlinlog(wlinlog)
        case('linw')     ! linear wall elements
          if (wlinlog.eq.2.or.wlinlog.eq.4) wlinlog=wlinlog-1
          call set_wlinlog(wlinlog)
        case('logxw')    ! log x wall elements plot
          if (wlinlog.le.2) wlinlog=wlinlog+2
          call set_wlinlog(wlinlog)
        case('linxw')    ! linear x wall elements plot
          if (wlinlog.gt.2) wlinlog=wlinlog-2
          call set_wlinlog(wlinlog)
        case('plwl')     ! plot in wall domain
          call walloff
          call plot_wall_elements
        case('sumw')     ! sum over walls
          stack_req=(/0,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          call sumw()
          graphlabel(iwwrk)=trim(graphlabel(iwwrk))//' summed '
          call set_wallbarlabel(wall_label(iwwrk))
        case('intw')     ! integrate over walls
          stack_req=(/0,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          call intw()
          graphlabel(iwwrk)=trim(graphlabel(iwwrk))//' integrated '
          call set_wallbarlabel(wall_label(iwwrk))
        case('intw*')    ! integrate over walls (no account of segments)
          stack_req=(/0,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          call intwstar()
          graphlabel(iwwrk)=trim(graphlabel(iwwrk))//' integrated '
          call set_wallbarlabel(wall_label(iwwrk))
        case('pswa')     ! per wall area
          stack_req=(/0,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          call pswa
        case('wclvs')    ! customise wall contour levels
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          stack_req=(/1,istack(iistack),0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call wclvs(istack(iistack),rstack,irstack)
          do iclv=1,istack(iistack)
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('wclvd')    ! customise wall contour level delta
          stack_req=(/0,3,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call wclvd(rstack(irstack-2),rstack(irstack-1),
     .               rstack(irstack))
          do iclv=1,3
            dummy=rpop(rstack,irstack)
          end do
        case('wclvn')    ! customise wall contour level range and number
          stack_req=(/1,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call wclvn(rstack(irstack-1),rstack(irstack),istack(iistack))
          do iclv=1,2
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('wclvp')    ! customise wall contour level by proportion
          stack_req=(/2,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          ix=ipop(istack,iistack)
          iy=ipop(istack,iistack)
          call wclvp(ix,iy)
          write(stmp,'(a,2x,a1,2i4,a2)')
     -     trim(wall_label(iwwrk)),'{',ix,iy,' }'
          wall_label(iwwrk)=stmp
          call set_wallbarlabel(wall_label(iwwrk))
        case('wwid')     ! wall line width factor
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call walllinewidth(rstack(irstack))
          wwidth = rstack(irstack)
          dummy=rpop(rstack,irstack)
        case('wabs')     ! abs() on wall stack
          stack_req=(/0,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do i=-1,nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=abs(wall_stack(i,is,iwwrk))
            end do
          end do
          wall_label(iwwrk)=trim(wall_label(iwwrk))//' (ABS)'
          call set_wallbarlabel(wall_label(iwwrk))
        case('ddw')      ! d/dw derivative along wall length
          stack_req=(/0,0,0,0,0,1,0/)
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
          wall_label(iwwrk)=trim(wall_label(iwwrk))//'/(m) (d/dw)'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wsumz')    ! sum over ns for wall stack data
          stack_req=(/2,0,0,0,0,1,0/)
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
     1         wall_stack(i,izstrt,iwwrk)+wall_stack(i,is,iwwrk)
            end do
          end do
          do i=1,nwdim(1,iwwrk)
            wall_stack(i,0,iwwrk)=
     1       wall_stack(i,izstrt,iwwrk)
          end do
          nwdim(2,iwwrk)=1
          if(izstrt .ne. izend) then
            wall_label(iwwrk)=
     1       trim(wall_label(iwwrk))//' SUM'//
     2       trim(textwwk(izstrt,iwwrk))//' to'//
     4       trim(textwwk(izend,iwwrk))
          else
            wall_label(iwwrk)=
     1       trim(wall_label(iwwrk))//trim(textwwk(izstrt,iwwrk))
          end if
          textwwk(0,iwwrk)=textbl(0)
          call pop_stack(stack_req)
          call set_wallbarlabel(wall_label(iwwrk))
        case('pwne')     ! wall data per ne
          stack_req=(/0,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          call pwne
        case('wm+')      ! wall stack addition
          stack_req=(/0,0,0,0,0,2,0/)
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
          stack_req=(/0,0,0,0,0,2,0/)
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
          stack_req=(/0,0,0,0,0,2,0/)
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
          stack_req=(/0,0,0,0,0,2,0/)
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
          stack_req=(/0,1,0,0,0,1,0/)
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
          stack_req=(/0,1,0,0,0,1,0/)
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
          stack_req=(/0,1,0,0,0,1,0/)
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
          stack_req=(/0,1,0,0,0,1,0/)
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
          stack_req=(/0,1,0,0,0,1,0/)
          if(bad_args(stack_req)) return   ! 5900
          do is=0,nwdim(2,iwwrk)-1
            do i=1,nwdim(1,iwwrk)
              if (rstack(irstack).ne.0.0_R8) then
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
          stack_req=(/0,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return   ! 5900
          iwwrk=iwwrk-1
        case('wdup')     ! duplicate top of wall stack
          stack_req=(/0,0,0,0,0,1,0/)
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
          stack_req=(/0,0,0,0,0,2,0/)
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
          stack_req=(/0,0,0,0,0,2,0/)
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
          if(iwwrk.eq.0) then
            write(*,*) 'No data on wall stack'
            return
          end if
          do i = iwwrk, 1, -1
            write(*,*) trim(wall_label(i)),
     &       (trim(textwwk(is,i)),is=0,nwdim(2,i)-1)
          end do
        case('wsel')     ! select wall element
          stack_req=(/1,0,0,0,0,1,0/)
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
          stack_req=(/1,0,0,0,0,1,0/)
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
          stack_req=(/0,0,0,0,0,2,0/)
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
          stack_req=(/0,0,0,0,0,2,0/)
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
          stack_req=(/0,1,0,0,0,1,0/)
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
        case('rwmma')    ! max(r1,w1)
          stack_req=(/0,1,0,0,0,1,0/)
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
        case('wmlt')     ! w<r
          stack_req=(/0,1,0,0,0,1,0/)
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
        case('wmle')     ! w<=r
          stack_req=(/0,1,0,0,0,1,0/)
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
        case('wmeq')     ! w=r
          stack_req=(/0,1,0,0,0,1,0/)
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
        case('wmne')     ! w!=r
          stack_req=(/0,1,0,0,0,1,0/)
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
        case('wmgt')     ! w>r
          stack_req=(/0,1,0,0,0,1,0/)
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
        case('wmge')     ! w=>r
          stack_req=(/0,1,0,0,0,1,0/)
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
          stack_req=(/0,0,0,0,0,1,0/)
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
          stack_req=(/0,0,0,0,0,1,0/)
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
           stack_req=(/0,0,0,0,0,1,0/)
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
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wexp')     ! exp(wall)
          stack_req=(/0,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=exp(wall_stack(i,is,iwwrk))
            end do
          end do
          wall_label(iwwrk)='Special'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wzrng')    ! select species range in wall stack
          stack_req=(/2,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          if(istack(iistack-1).gt.istack(iistack)) then
            write(*,*) 'Parameters for wzrng must be given in ',
     &       'ascending order!'
            return
          end if
          if(0.gt.istack(iistack-1).or.
     &       istack(iistack-1).ge.nwdim(2,iwwrk)) then
            write(*,*) 'Faulty first parameter is = ',istack(iistack-1),
     &       ' for wzrng!'
            return
          endif
          if(0.gt.istack(iistack).or.
     &       istack(iistack).ge.nwdim(2,iwwrk)) then
            write(*,*) 'Faulty second parameter is = ',istack(iistack),
     &       ' for wzrng!'
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
          stack_req=(/0,0,0,0,0,1,0/)
          if(bad_args(stack_req)) return
          do is=0,nwdim(2,iwwrk)-1
            do i = 1, nwdim(1,iwwrk)
              wall_stack(i,is,iwwrk)=-wall_stack(i,is,iwwrk)
            end do
          end do
          wall_label(iwwrk)=trim(wall_label(iwwrk))//' (CHS)'
          call set_wallbarlabel(wall_label(iwwrk))
        case('wtos')     ! wall array -> m
          stack_req=(/0,0,0,0,0,1,0/)
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
        case('wtot')     ! wall array --> target
          stack_req=(/0,0,0,0,0,1,0/)
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
          stack_req=(/0,0,0,0,0,1,0/)
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
          stack_req=(/1,0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          arg=spop(sstack,isstack)
          ndata=ipop(istack,iistack)
          nz=1
          open(99,file=arg,status='old')
          do i=1,ndata
            read(99,*) r1(i),r2(i),z1(i),z2(i),data(i)
          enddo
          values(5)=min(minval(r1),minval(r2))                   !rmin
          values(6)=max(maxval(r1),maxval(r2))                   !rmax
          values(7)=min(minval(z1),minval(z2))                   !zmin
          values(8)=max(maxval(z1),maxval(z2))                   !zmax
          write(*,*) 'values',(i,values(i),i=5,8)
          call wall_copy(data,1.0_R8,nz,textbl)
          wall_label(iwwrk)=arg
          call set_wallbarlabel(wall_label(iwwrk))
        case('wallf')    ! add experimental wall data to plot
          stack_req=(/0,0,1,0,0,0,0/)
          if(bad_args(stack_req)) return
          wallfile=spop(sstack,isstack)
          call plot_wallfile
        case('outpw')    ! prints the value of a specific wall element
          stack_req=(/1,0,1,0,0,1,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          i=ipop(istack,iistack)
          do is=0,nwdim(2,iwwrk)-1
            if(nwdim(2,iwwrk) .eq. 1) then
              write(*,*) trim(stmp),wall_stack(i,is,iwwrk)
              if(writing_to_file) then
                write(nwtf,*) trim(stmp),wall_stack(i,is,iwwrk)
              endif
            else
              write(*,*) trim(stmp),is,wall_stack(i,is,iwwrk)
              if(writing_to_file) then
                write(nwtf,*) trim(stmp),is,wall_stack(i,is,iwwrk)
              endif
            endif
          enddo
        case('wprt')     ! output the top of wall stack to a file
          stack_req=(/0,0,1,0,0,1,0/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          open(99,file=stmp)
          do is=0,nwdim(2,iwwrk)-1
            write(99,'(''#'',1x,a,1x,a)')
     1       trim(wall_label(iwwrk)),trim(textwwk(is,iwwrk))
            do i=1,nbnd
              write(99,*) i,is,wall_stack(i,is,iwwrk)
            enddo
          enddo
          close(99)
        case('wtab')     ! output a table from the top of wall stack to a file
          stack_req=(/2,0,1,0,0,1,0/)
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
          write(hlp_frm,'(a,i3,a)') '(a1,1x,a,1x,',nz,'a13)'
          write(*,*) hlp_frm
          write(99,hlp_frm) '#',
     1       trim(wall_label(iwwrk)),(trim(textwwk(is,iwwrk)),is=0,nz-1)
          write(hlp_frm,'(a,i3,a)') '(i4,1p,',nz,'(1e12.5,1x))'
          write(*,*) hlp_frm
          do ix=ix1,ix2
            write(99,hlp_frm) ix,wall_stack(ix,0:nz-1,iwwrk)
          enddo
          close(99)
        case('wlld')     ! ANK wall loading
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call b2pwld(istack(iistack),BoRiS)
          call pop_stack(stack_req)
        case('woff')     ! turn off wall plots
          call walloff
        case default WFVALUE
          cflag=.false.
      END SELECT WFVALUE
c
      goto 5900

 5059 continue
      stack_req=(/0,0,0,0,0,1,0/)
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
      integer i,ix,iy,ibnd,iclv,is,iseg

      cflag=.true.
      TFVALUE: select case (cval)
        case('logtrg')   ! log target elements
          if (mod(trglinlog,2).eq.1) trglinlog=trglinlog+1
          call set_trglinlog(trglinlog)
        case('lintrg')   ! linear target elements
          if (trglinlog.eq.2.or.trglinlog.eq.4) trglinlog=trglinlog-1
          call set_trglinlog(trglinlog)
        case('logxtrg')  ! log x target elements plot
          if (trglinlog.le.2) trglinlog=trglinlog+2
          call set_trglinlog(trglinlog)
        case('linxtrg')  ! linear x target elements plot
          if (trglinlog.gt.2) trglinlog=trglinlog-2
          call set_trglinlog(trglinlog)
        case('tclvs')    ! customize target contour levels
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          stack_req=(/1,istack(iistack),0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call tclvs(istack(iistack),rstack,irstack)
          do iclv=1,istack(iistack)
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('tclvd')    ! customise target contour level delta
          stack_req=(/0,3,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call tclvd
     .     (rstack(irstack-2),rstack(irstack-1),rstack(irstack))
          do iclv=1,3
            dummy=rpop(rstack,irstack)
          end do
        case('tclvn')    ! customise target contour level range and number
          stack_req=(/1,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call tclvn(rstack(irstack-1),rstack(irstack),istack(iistack))
          do iclv=1,2
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('tclvp')    ! customise target contour level by proportion
          stack_req=(/2,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          ix=ipop(istack,iistack)
          iy=ipop(istack,iistack)
          call tclvp(ix,iy)
          write(stmp,'(a,2x,a1,2i4,a2)')
     -     trim(target_label(itwrk)),'{',ix,iy,' }'
          target_label(itwrk)=stmp
          call set_targetbarlabel(target_label(itwrk))
        case('trg.w')    ! line plot versus target wall
          stack_req=(/1,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          if(1.le.istack(iistack).and.istack(iistack).le.ndepth_nml)
     >     then
            call trgw(istack(iistack))
          else
            write(*,*)
     .       'Faulty iz = ',istack(iistack),' chosen for trg.w!'
          endif
          call pop_stack(stack_req)
        case('trg.z')    ! line plot versus target depth
          stack_req=(/1,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          if(1.le.istack(iistack).and.istack(iistack).le.nbnd) then
            call trgz(istack(iistack))
          else
            write(*,*)
     .       'Faulty ibnd = ',istack(iistack),' chosen for trg.z!'
          endif
          call pop_stack(stack_req)
        case('sumtrgw')  ! sum over targets along walls
          stack_req=(/0,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call sumtrgw()
          target_label(itwrk)=
     1     trim(target_label(itwrk))//' summed over w'
          call set_targetbarlabel(target_label(itwrk))
        case('sumtrgz')  ! sum over targets along z
          stack_req=(/0,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call sumtrgz()
          target_label(itwrk)=
     1     trim(target_label(itwrk))//' summed over z'
          call set_targetbarlabel(target_label(itwrk))
        case('inttrgw')  ! integrate over targets along walls
          stack_req=(/0,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call inttrgw()
          target_label(itwrk)=
     1     trim(target_label(itwrk))//' integrated along w'
          call set_targetbarlabel(target_label(itwrk))
        case('inttrgw*') ! integrate over targets along walls (no account of segments)
          stack_req=(/0,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call inttrgwstar()
          target_label(itwrk)=
     1     trim(target_label(itwrk))//' integrated along w'
          call set_targetbarlabel(target_label(itwrk))
        case('inttrgz')  ! integrate over targets along z
          stack_req=(/0,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call inttrgz()
          target_label(itwrk)=
     1     trim(target_label(itwrk))//' integrated along z'
          call set_targetbarlabel(target_label(itwrk))
        case('trgabs')   ! absolute value of target data
          stack_req=(/0,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          target_stack(1:nbnd,:,itwrk)=abs(target_stack(1:nbnd,:,itwrk))
          target_label(itwrk)=trim(target_label(itwrk))//' (ABS) '
          call set_targetbarlabel(target_label(itwrk))
        case('trgddw')   ! target d/dw derivative along wall length
          stack_req=(/0,0,0,0,0,0,1/)
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
          target_label(itwrk)=trim(target_label(itwrk))//'(/m) (d/dw)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgddp')   ! target d/dp derivative along perpendicular direction
          stack_req=(/0,0,0,0,0,0,1/)
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
          target_label(itwrk)=trim(target_label(itwrk))//'(/m) (d/dp)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgddz')   ! target d/dz derivative along target depth
          stack_req=(/0,0,0,0,0,0,1/)
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
          target_label(itwrk)=trim(target_label(itwrk))//'(/m) (d/dz)'
          call set_targetbarlabel(target_label(itwrk))
        case('rtrgm+')   ! target stack + real
          stack_req=(/0,1,0,0,0,0,1/)
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
        case('rtrgm-')   ! target stack - real
          stack_req=(/0,1,0,0,0,0,1/)
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
        case('rtrgm*')   ! target stack * real
          stack_req=(/0,1,0,0,0,0,1/)
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
        case('rtrgm/')   ! target stack / real
          stack_req=(/0,1,0,0,0,0,1/)
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
        case('rtrgm**')  ! target stack ** real
          stack_req=(/0,1,0,0,0,0,1/)
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
        case('trgtemp')  ! target temperature (K)
          call target_copy(target_temp,1.0_R8)
          target_label(itwrk)='Target temperature (K)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgdrop')  ! drop top of target stack
          stack_req=(/0,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          itwrk=itwrk-1
        case('trgwsel')  ! select wall element in target stack
          stack_req=(/1,0,0,0,0,0,1/)
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
        case('trgzsel')  ! select depth element in target stack
          stack_req=(/1,0,0,0,0,0,1/)
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
        case('rtrgmmi')  ! min(r1,trg)
          stack_req=(/0,1,0,0,0,0,1/)
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
        case('rtrgmma')  ! max (r1,trg)
          stack_req=(/0,1,0,0,0,0,1/)
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
          stack_req=(/0,1,0,0,0,0,1/)
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
          stack_req=(/0,1,0,0,0,0,1/)
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
          stack_req=(/0,1,0,0,0,0,1/)
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
          stack_req=(/0,1,0,0,0,0,1/)
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
          stack_req=(/0,1,0,0,0,0,1/)
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
          stack_req=(/0,1,0,0,0,0,1/)
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
        case('trglog')   ! log10(target)
          stack_req=(/0,0,0,0,0,0,1/)
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
        case('trgln')    ! ln(target)
          stack_req=(/0,0,0,0,0,0,1/)
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
        case('trgsqrt')  ! sqrt(target)
          stack_req=(/0,0,0,0,0,0,1/)
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
        case('trgexp')   ! exp(target)
          stack_req=(/0,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              target_stack(ibnd,i,itwrk)=
     1         exp(target_stack(ibnd,i,itwrk))
            end do
          end do
          target_label(itwrk)='Special'
          call set_targetbarlabel(target_label(itwrk))
        case('trgchs')   ! change sign in target stack
          stack_req=(/0,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          do i=1,ndepth_nml
            do ibnd = 1, nbnd
              target_stack(ibnd,i,itwrk)=-target_stack(ibnd,i,itwrk)
            end do
          end do
          target_label(itwrk)=trim(target_label(itwrk))//' (CHS)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgtow')   ! target array --> wall
          stack_req=(/0,0,0,0,0,0,1/)
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
        case('trgpos')   ! target stack with 1 if element > 0, else 0
          stack_req=(/0,0,0,0,0,0,1/)
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
        case('outptrg')  ! prints the value of a specific target element
          stack_req=(/2,0,1,0,0,0,1/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          i=ipop(istack,iistack)
          ibnd=ipop(istack,iistack)
          write(*,*) trim(stmp),target_stack(ibnd,i,itwrk)
          if(writing_to_file) then
            write(nwtf,*) trim(stmp),target_stack(ibnd,i,itwrk)
          endif
        case('trgprt')   ! output the target stack to a file
          stack_req=(/0,0,1,0,0,0,1/)
          if(bad_args(stack_req)) return
          stmp=spop(sstack,isstack)
          open(99,file=stmp)
          write(99,'(''#'',1x,a)') trim(target_label(itwrk))
          do i=1,ndepth_nml
            do ibnd=1,nbnd
              write(99,*) ibnd,i,target_stack(ibnd,i,itwrk)
            enddo
          enddo
          close(99)
        case('trgtab')   ! output a table from the top of target stack to a file
          stack_req=(/4,0,1,0,0,0,1/)
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
          write(99,*) '#', trim(target_label(itwrk))
          do iy=iy1,iy2
            do ix=ix1,ix2
              write(99,'(2i4,1pe12.5)') ix,iy,target_stack(ix,iy,itwrk)
            enddo
          enddo
          close(99)
        case('trgsurf')  ! plot target plot
          stack_req=(/0,0,0,0,0,0,1/)
          if(bad_args(stack_req)) return
          call trgsurf(1)
        case('kappa_ab') ! target heat conductivity computed as 1/(AT+B)
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              if (target_depth(ibnd,i).ge.coating_thickness(ibnd)) then
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
        case('kappa')    ! target heat conductivity computed from polynomial fit
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              trg_tmp(ibnd,i) = Kappa(ibnd,i,target_temp(ibnd,i))
            end do
          end do
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='Target heat conductivity (W/m/K)'
          call set_targetbarlabel(target_label(itwrk))
        case('dkappadt') ! target heat conductivity derivative
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
        case('kappapdt') ! perp. target heat conductivity derivative
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
        case('cp')       ! target specific heat capacity
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              trg_tmp(ibnd,i)=Cp(ibnd,i,target_temp(ibnd,i))
            end do
          end do
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='Target heat capacity (J/kg/K)'
          call set_targetbarlabel(target_label(itwrk))
        case('rho')      ! target mass density
          do ibnd = 1, nbnd
            do i = 1, ndepth_nml
              trg_tmp(ibnd,i)=mass_density(ibnd,i)
            end do
          end do
          call target_copy(trg_tmp,1.0_R8)
          call customise_label
     .     ('Target mass density (kg/m^{3})',target_label(itwrk))
          call set_targetbarlabel(target_label(itwrk))
        case('trgwidth') ! target element width
          call target_copy(target_width,1.0_R8)
          target_label(itwrk)='Target element width (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgdepth') ! target element depth position
          call target_copy(target_depth,1.0_R8)
          target_label(itwrk)='Target element depth (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgvol')   ! target element volume
          call target_copy(target_vol,1.0_R8)
          call customise_label
     .     ('Target element volume (m^{3})',target_label(itwrk))
          call set_targetbarlabel(target_label(itwrk))
        case('trgdup')   ! duplicate top of target stack
          stack_req=(/0,0,0,0,0,0,1/)
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
        case('trgswap')  ! swap top two elements of target stack
          stack_req=(/0,0,0,0,0,0,2/)
          if(bad_args(stack_req)) return
          do ibnd=1,nbnd
            do i=1,ndepth_nml
              call rswap(target_stack(ibnd,i,itwrk),
     1                   target_stack(ibnd,i,itwrk-1))
            end do
          end do
          call strswap(target_label(itwrk),target_label(itwrk-1))
          call set_targetbarlabel(target_label(itwrk))
        case('trgover')  ! t2 t1 --> t2 t1 t2
          stack_req=(/0,0,0,0,0,0,2/)
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
        case('trglist')  ! list contents of target stack
          if(itwrk.eq.0) then
            write(*,*) 'No data on target stack'
            return
          end if
          do i = itwrk, 1, -1
            write(*,'(a)') trim(target_label(i))
          end do
        case('trgmmin')  ! min(t1,t2)
          stack_req=(/0,0,0,0,0,0,2/)
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
        case('trgmmax')  ! max(t1,t2)
          stack_req=(/0,0,0,0,0,0,2/)
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
        case('trgr1')    ! target element R1 position (m)
          call target_copy(targetcr(:,:,1),1.0_R8)
          target_label(itwrk)='Target element R1 position (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgr2')    ! target element R2 position (m)
          call target_copy(targetcr(:,:,2),1.0_R8)
          target_label(itwrk)='Target element R2 position (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgz1')    ! target element Z1 position (m)
          call target_copy(targetcz(:,:,1),1.0_R8)
          target_label(itwrk)='Target element Z1 position (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgz2')    ! target element Z2 position (m)
          call target_copy(targetcz(:,:,2),1.0_R8)
          target_label(itwrk)='Target element Z2 position (m)'
          call set_targetbarlabel(target_label(itwrk))
        case('trgnorth') ! North boundary
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
        case('trgsouth') ! South boundary
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
        case('trgseg')   ! target segment
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('trgibnd')  ! trgibnd
          do ibnd = 1, nbnd
            trg_tmp(ibnd,1:ndepth_nml) = real(ibnd)
          end do
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='Target element index'
          call set_targetbarlabel(target_label(itwrk))
        case ('trgiz')    ! trgiz
          do is = 1, ndepth_nml
            trg_tmp(1:nbnd,is) = real(is)
          end do
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='Target depth index'
          call set_targetbarlabel(target_label(itwrk))
        case ('itrglt')   ! itrglt
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('itrgle')   ! itrgle
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          is = istack(iistack)
          if(is.gt.nbnd .or. is.le.0) then
            write(*,*) 'No such wall element'
            return
          endif
          trg_tmp(1:is,1:ndepth_nml) = 1.0_R8
          trg_tmp(is+1:nbnd,1:ndepth_nml) = 0.0_R8
          write(*,*) '"itrgle" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='itrgle'
          call set_targetbarlabel(target_label(itwrk))
        case ('itrgeq')   ! itrgeq
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('itrgne')   ! itrgne
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          is = istack(iistack)
          if(is.gt.nbnd .or. is.le.0) then
            write(*,*) 'No such wall element'
            return
          endif
          trg_tmp(1:nbnd,1:ndepth_nml) = 1.0_R8
          trg_tmp(is,1:ndepth_nml) = 0.0_R8
          write(*,*) '"itrgne" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='itrgeq'
          call set_targetbarlabel(target_label(itwrk))
        case ('itrggt')   ! itrggt
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('itrgge')   ! itrgge
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          is = istack(iistack)
          if(is.gt.nbnd .or. is.le.0) then
            write(*,*) 'No such wall element'
            return
          endif
          trg_tmp(1:is-1,1:ndepth_nml) = 0.0_R8
          trg_tmp(is:nbnd,1:ndepth_nml) = 1.0_R8
          write(*,*) '"itrgge" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='itrgge'
          call set_targetbarlabel(target_label(itwrk))
        case ('ztrglt')   ! ztrglt
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('ztrgle')   ! ztrgle
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp = istack(iistack)
          if(itmp.gt.ndepth_nml .or. itmp.le.0) then
            write(*,*) 'No such wall depth index'
            return
          endif
          trg_tmp(1:nbnd,1:itmp) = 1.0_R8
          trg_tmp(1:nbnd,itmp+1:ndepth_nml) = 0.0_R8
          write(*,*) '"ztrgle" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='ztrgle'
          call set_targetbarlabel(target_label(itwrk))
        case ('ztrgeq')   ! ztrgeq
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('ztrgne')   ! ztrgne
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp = istack(iistack)
          if(itmp.gt.ndepth_nml .or. itmp.le.0) then
            write(*,*) 'No such wall depth index'
            return
          endif
          trg_tmp(1:nbnd,1:ndepth_nml) = 1.0_R8
          trg_tmp(1:nbnd,itmp) = 0.0_R8
          write(*,*) '"ztrgne" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='ztrgeq'
          call set_targetbarlabel(target_label(itwrk))
        case ('ztrggt')   ! ztrggt
          stack_req=(/1,0,0,0,0,0,0/)
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
        case ('ztrgge')   ! ztrgge
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          itmp = istack(iistack)
          if(itmp.gt.ndepth_nml .or. itmp.le.0) then
            write(*,*) 'No such wall depth index'
            return
          endif
          trg_tmp(1:nbnd,1:itmp-1) = 0.0_R8
          trg_tmp(1:nbnd,itmp:ndepth_nml) = 1.0_R8
          write(*,*) '"ztrgge" option with ',ipop(istack,iistack)
          call target_copy(trg_tmp,1.0_R8)
          target_label(itwrk)='ztrgge'
          call set_targetbarlabel(target_label(itwrk))
        case('trgm+')    ! target stack addition
          stack_req=(/0,0,0,0,0,0,2/)
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
        case('trgm-')    ! target stack subtraction
          stack_req=(/0,0,0,0,0,0,2/)
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
        case('trgm*')    ! target stack multiplication
          stack_req=(/0,0,0,0,0,0,2/)
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
        case('trgm/')    ! target stack division
          stack_req=(/0,0,0,0,0,0,2/)
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
        case('me')       ! electron mass
          stack_name='real'
          call rpush(rstack,irstack,me)
        case('mp')       ! proton mass
          stack_name='real'
          call rpush(rstack,irstack,mp)
        case('ev')       ! electron charge
          stack_name='real'
          call rpush(rstack,irstack,ev)
        case('qe')       ! electron charge
          stack_name='real'
          call rpush(rstack,irstack,qe)
        case('mu0')      ! magnetic permittivity of the vacuum
          stack_name='real'
          call rpush(rstack,irstack,mu0)
        case('eps0')     ! electric permeability of the vacuum
          stack_name='real'
          call rpush(rstack,irstack,eps0)
        case('pi')       ! pi
          stack_name='real'
          call rpush(rstack,irstack,pi)
        case('c')        ! speed of light
          stack_name='real'
          call rpush(rstack,irstack,c)
        case('avogr')    ! Avogadro's number
          stack_name='real'
          call rpush(rstack,irstack,avogr)
        case('kbolt')    ! Boltzmann constant
          stack_name='real'
          call rpush(rstack,irstack,kbolt)
        case('jxi')      ! jxi
          stack_name='integer'
          call ipush(istack,iistack,jxi)
        case('jxa')      ! jxa
          stack_name='integer'
          call ipush(istack,iistack,jxa)
        case('jsep')     ! jsep
          stack_name='integer'
          call ipush(istack,iistack,jsep)
        case('ixref')    ! ixref
          stack_name='integer'
          call ipush(istack,iistack,ixref)
        case('iyref')    ! iyref
          stack_name='integer'
          call ipush(istack,iistack,iyref)
        case('nx')       ! nx
          stack_name='integer'
          call ipush(istack,iistack,nx)
        case('ny')       ! ny
          stack_name='integer'
          call ipush(istack,iistack,ny)
        case('ns')       ! ns
          stack_name='integer'
          call ipush(istack,iistack,ns)
        case('nspecies') ! nspecies
          stack_name='integer'
          call ipush(istack,iistack,nspecies)
        case('natm')     ! natm
          stack_name='integer'
          call ipush(istack,iistack,natm)
        case('nmol')     ! nmol
          stack_name='integer'
          call ipush(istack,iistack,nmol)
        case('nion')     ! nion
          stack_name='integer'
          call ipush(istack,iistack,nion)
        case('nxtl')     ! nxtl
          stack_name='integer'
          if (nnreg(0).ge.7) then
            call ipush(istack,iistack,nxtl)
          else
            write(*,*) cval,
     .       'option available for DN cases only: Nothing done.'
          endif
        case('nxtr')     ! nxtr
          stack_name='integer'
          if (nnreg(0).ge.7) then
            call ipush(istack,iistack,nxtr)
          else
            write(*,*) cval,
     .       'option available for DN cases only: Nothing done.'
          endif
        case('nbnd')     ! nbnd
          stack_name='integer'
          call ipush(istack,iistack,nbnd)
        case('ndepth')   ! ndepth_nml
          stack_name='integer'
          call ipush(istack,iistack,ndepth_nml)
        case('ismain')   ! ismain
          stack_name='integer'
          call ipush(istack,iistack,ismain)
        case('ismain0')  ! ismain0
          stack_name='integer'
          call ipush(istack,iistack,ismain0)
        case('ns_ext')   ! ns_ext
          stack_name='integer'
          call ipush(istack,iistack,ns_ext)
        case default CNVALUE
          cflag=.false.
      END SELECT CNVALUE

      if(cflag) write(*,*) 'Copied constant ',trim(cval),' to ',
     . trim(stack_name),' stack.'
      return

      end subroutine select_constants

      subroutine select_flvalue
      integer is, iclv

      cflag=.true.
      FLVALUE: select case(cval)
        case('f30')      ! write fort.30
          call write_f30(nx, ny)
        case('f31')      ! write fort.31
          call write_f31(nx, ny)
        case('glob')     ! transform an extra label into a global header
          call set_globalheader
        case('both')     ! both phys and comp
          typeofplot=3   ! both
        case('phys')     ! physical domain only
          typeofplot=1   ! physical domain only
        case('comp')     ! computational domain only
          typeofplot=2   ! computational domain only
        case('trigmesh') ! overlay triangular grid
          call trigmesh
        case('mesh')     ! overlay mesh
          call mesh
        case('expl')     ! toggle value of expline
          expline=.not.expline
        case('page')     ! set number of graphs per page
          stack_req=(/2,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          iy=ipop(istack,iistack)
          ix=ipop(istack,iistack)
          if(ix.le.0.or.iy.le.0) then
            write(*,*) 'Invalid values for page!'
            write(*,*) 'Must be strictly positive'
            return
          end if
          call b2page(ix,iy)
        case('logf')     ! log f plot
          if (mod(linlog,2).eq.1) linlog=linlog+1
        case('linf')     ! linear f plot
          if (linlog.eq.2.or.linlog.eq.4) linlog=linlog-1
        case('logx')     ! log x plot
          if(linlog.le.2) linlog=linlog+2
        case('linx')     ! linear x plot
          if(linlog.gt.2) linlog=linlog-2
        case('clvs')     ! customise contour levels
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          stack_req=(/1,istack(iistack),0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call clvs(istack(iistack),rstack,irstack)
          do iclv=1,istack(iistack)
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('clvd')     ! customise contour level delta
          stack_req=(/0,3,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call clvd(rstack(irstack-2),rstack(irstack-1),rstack(irstack))
          do iclv=1,3
            dummy=rpop(rstack,irstack)
          end do
        case('clvn')     ! customise contour level range and number
          stack_req=(/1,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call clvn(rstack(irstack-1),rstack(irstack),istack(iistack))
          do iclv=1,2
            dummy=rpop(rstack,irstack)
          end do
          call pop_stack(stack_req)
        case('clvp')     ! customise contour level by proportion
          stack_req=(/2,0,0,1,0,0,0/)
          if(bad_args(stack_req)) return
          ix=ipop(istack,iistack)
          iy=ipop(istack,iistack)
          call clvp(ix,iy)
          write(stmp,'(a,2x,a1,2i4,a2)')
     -     trim(graphlabel(iwrk)),'{',ix,iy,' }'
          graphlabel(iwrk)=stmp
        case('debu','debug')  ! set debugging level
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          b2plot_debug=ipop(istack,iistack)
          write(*,*) 'Debugging set to ',b2plot_debug
        case('exps')     ! set symbol for experimental data
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          expsymbol=ipop(istack,iistack)
        case('eclv')     ! do not fill contour Levels
          call nofill(.true.)
        case('fclv')     ! do fill contour Levels
          call nofill(.false.)
        case('font')     ! set font
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call font(istack(iistack))
          call pop_stack(stack_req)
        case('lacl')     ! label orientation along contour levels
          call lblalongcl(.true.)
        case('lhor')     ! label orientation horizontal
          call lblalongcl(.false.)
        case('lbon')     ! labeling on / off
          call labelon(.true.)
        case('lbof')     ! labeling on / off
          call labelon(.false.)
        case('lwid')     ! line width factor
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call linewidth(rstack(irstack))
          lwidth=rpop(rstack,irstack)
        case('zwid')     ! axis line width factor
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call axislinewidth(rstack(irstack))
          zwidth=rpop(rstack,irstack)
        case('lfsz')     ! line label font size
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_line_font_size(istack(iistack))
          lfsz=istack(iistack)
          call pop_stack(stack_req)
        case('grid')     ! grid on
          call gridon(.true.)
        case('trig')     ! triangular grid on
          call trigon(.true.)
        case('clvl')     ! contour level lines off
          call clvl(.false.)
        case('spln')     ! switch on interpolation
          call ipolate(.true.)
        case('arow')     ! change default arrow length
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          arrowl=rpop(rstack,irstack)
        case('earow')    ! change default arrow length for Eirene stack data
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          earrowl=rpop(rstack,irstack)
        case('arww')     ! change default arrow width
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          arroww=rpop(rstack,irstack)
        case('acol')     ! colour index for arrows
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          acol=istack(iistack)
          call pop_stack(stack_req)
        case('sep')      ! plot separatrix
          plotting_sep=.not.plotting_sep
          call separatrixon(.true.)
        case('vesl')     ! plot vessel walls from mesh.extra file
          call vesselon(.true.)
        case('mmid')     ! map to midplane
          call set_midplane_mapping(.true.)
        case('fpsi')     ! plot versus full psi
          call set_full_psi_plotting(.true.)
        case('mpsi')     ! plot versus normalised psi
          call set_norm_psi_plotting(.true.)
        case('rpsi')     ! plot versus root(normalised psi)
          call set_root_psi_plotting(.true.)
        case('nocp')     ! do not plot contour plot
          call contourplot(.false.)
        case('nice')     ! toggle whether axes will have nice values
          call adjust_axes(.true.)
        case('bold')     ! toggle whether to bold the axis labels
          call adjust_bold(.true.)
        case('sepw')     ! change separatrix line width
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          sepwidth=rpop(rstack,irstack)
        case('ewid')     ! change experimental data line width
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          ewidth=rpop(rstack,irstack)
        case('owid')     ! change outline line width
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          owidth=rpop(rstack,irstack)
        case('vwid')     ! change vessel line width
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          vwidth=rpop(rstack,irstack)
        case('sq')       ! square aspect ratio
          call set_square
        case('xpol')     ! use poloidal distance for f.x
         xparallel=.false.
        case('xpar')     ! use parallel distance for f.x
         xparallel=.true.
        case('a4p')      ! A4 portrait
          call set_a4p
        case('a4l')      ! A4 landscape
          call set_a4l
        case('same')     ! next graph on same axes
          call set_same_page
        case('filo')     ! fill option
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_fill_option(istack(iistack))
          call pop_stack(stack_req)
        case('i+')       ! i1 + i2
          stack_req=(/2,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          it2=ipop(istack,iistack)
          it1=ipop(istack,iistack)
          call ipush(istack,iistack,it1+it2)
        case('i-')       ! i1 - i2
          stack_req=(/2,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          it2=ipop(istack,iistack)
          it1=ipop(istack,iistack)
          call ipush(istack,iistack,it1-it2)
        case('i*')       ! i1 * i2
          stack_req=(/2,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          it2=ipop(istack,iistack)
          it1=ipop(istack,iistack)
          call ipush(istack,iistack,it1*it2)
        case('i/')       ! i1 / i2
          stack_req=(/2,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          it2=ipop(istack,iistack)
          it1=ipop(istack,iistack)
          if (it2.eq.0) then
            write(*,*) 'Division by zero : Operation not permitted'
            return
          endif
          call ipush(istack,iistack,it1/it2)
        case('i**')      ! i1 ** i2
          stack_req=(/2,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          it2=ipop(istack,iistack)
          it1=ipop(istack,iistack)
          if (it1.eq.0 .and. it2.eq.0) then
            write(*,*) 'Zero to the zeroth power : Undefined'
            return
          endif
          call ipush(istack,iistack,it1**it2)
        case('itor')     ! integer --> real stack
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          dummy=real(ipop(istack,iistack),R8)
          call rpush(rstack,irstack,dummy)
        case('itostr')   ! integer --> string stack
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          write(hlp_frm,'(a,i2,a)')
     .     '(i',len_of_digits(istack(iistack)),')'
          write(stmp,hlp_frm) ipop(istack,iistack)
          call spush(sstack,isstack,stmp)
        case('itom')     ! integer --> matrix stack
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(iwrk .eq. wklng) call dropstackbottom()
          iwrk=iwrk+1
          ndim(1,iwrk)=nx
          ndim(2,iwrk)=ny
          ndim(3,iwrk)=1
          work(-1:nx,-1:ny,0,iwrk)=real(istack(iistack),R8)
          values(9)=real(istack(iistack),R8)
          values(10)=real(istack(iistack),R8)
          textwk(0,iwrk)=textbl(0)
          write(hlp_frm,'(a,i2,a)')
     .     '(i',len_of_digits(istack(iistack)),')'
          write(graphlabel(iwrk),hlp_frm) ipop(istack,iistack)
        case('itoe')     ! integer --> Eirene stack
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
#ifdef B25_EIRENE
          if(iewrk .eq. wklng) call dropeirenestackbottom()
          iewrk=iewrk+1
          nedim(1,iewrk)=ntrii
          nedim(2,iewrk)=1
          eirene_stack(1:ntrii,0,iewrk)=real(istack(iistack),R8)
          values(38)=real(istack(iistack),R8)
          values(39)=real(istack(iistack),R8)
          textewk(0,iewrk)=textbl(0)
          write(hlp_frm,'(a,i2,a)')
     .     '(i',len_of_digits(istack(iistack)),')'
          write(eirene_label(iewrk),hlp_frm) ipop(istack,iistack)
#else
          write(*,*) 'Must compile with Eirene for itoe function'
          call pop_stack(stack_req)
#endif
        case('itow')     ! integer --> wall stack
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(iwwrk .eq. wklng) call dropwallstackbottom()
          iwwrk=iwwrk+1
          nwdim(1,iwwrk)=nbnd
          nwdim(2,iwwrk)=1
          wall_stack(1:nbnd,0,iwwrk)=real(istack(iistack),R8)
          values(14)=real(istack(iistack),R8)
          values(15)=real(istack(iistack),R8)
          textwwk(0,iwwrk)=textbl(0)
          write(hlp_frm,'(a,i2,a)')
     .     '(i',len_of_digits(istack(iistack)),')'
          write(wall_label(iwwrk),hlp_frm) ipop(istack,iistack)
          call set_wallbarlabel(wall_label(iwwrk))
        case('itot')     ! integer --> target stack
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(itwrk .eq. wklng) call droptargetstackbottom()
          itwrk=itwrk+1
          target_stack(1:nbnd,1:ndepth_nml,itwrk)=
     .               real(istack(iistack),R8)
          values(32)=real(istack(iistack),R8)
          values(33)=real(istack(iistack),R8)
          write(hlp_frm,'(a,i2,a)')
     .     '(i',len_of_digits(istack(iistack)),')'
          write(target_label(itwrk),hlp_frm) ipop(istack,iistack)
          call set_targetbarlabel(target_label(itwrk))
        case('iprt')     ! output the integer stack
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is = iistack, 1, -1
            write(*,*) istack(is)
            if(writing_to_file) then
              write(nwtf,*) istack(is)
            endif
          enddo
        case('r+')       ! r1 + r2
          stack_req=(/0,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          t2=rpop(rstack,irstack)
          t1=rpop(rstack,irstack)
          call rpush(rstack,irstack,t1+t2)
        case('r-')       ! r1 - r2
          stack_req=(/0,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          t2=rpop(rstack,irstack)
          t1=rpop(rstack,irstack)
          call rpush(rstack,irstack,t1-t2)
        case('r*')       ! r1 * r2
          stack_req=(/0,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          t2=rpop(rstack,irstack)
          t1=rpop(rstack,irstack)
          call rpush(rstack,irstack,t1*t2)
        case('r/')       ! r1 / r2
          stack_req=(/0,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          t2=rpop(rstack,irstack)
          t1=rpop(rstack,irstack)
          if (t2.eq.0.0_R8) then
            write(*,*) 'Division by zero : Operation not permitted'
            return
          endif
          call rpush(rstack,irstack,t1/t2)
        case('r**')      ! r1 ** r2
          stack_req=(/0,2,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          t2=rpop(rstack,irstack)
          t1=rpop(rstack,irstack)
          if (t1.eq.0.0_R8 .and. t2.eq.0.0_R8) then
            write(*,*) 'Zero to the zeroth power : Undefined'
            return
          endif
          call rpush(rstack,irstack,t1**t2)
        case('rtostr')   ! real --> string stack
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          write(stmp,'(1pe15.8)') rpop(rstack,irstack)
          call spush(sstack,isstack,stmp)
        case('rprt')     ! output the real stack
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          do is = irstack, 1, -1
            write(*,'(1pe15.8)') rstack(is)
            if(writing_to_file) then
              write(nwtf,'(1pe15.8)') rstack(is)
            endif
          enddo
        case('new')      ! new page
          call set_new_page
        case('lbsz')     ! label size
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_lbsz(rpop(rstack,irstack))
        case('glsz')     ! globalheader size
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_glsz(rpop(rstack,irstack))
        case('ndec')     ! number of decades
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(istack(iistack).eq.0) then  ! restore default
            call set_ndec(10)
          else if(istack(iistack).gt.0) then
            call set_ndec(istack(iistack))
          else
            write(6,*) 'A positive number of decades, please!'
          endif
          call pop_stack(stack_req)
        case('aspc')     ! preserve aspect ratio flag
          call swap_aspect
        case('col')      ! colour
#ifndef PLPLOT
          call sfseti('TY',0)
#endif
          color_set=.true.
        case('bw')       ! black and white
#ifndef PLPLOT
          call sfseti('TY',-4)
#endif
          color_set=.false.
        case('logt')     ! log time for trajectories
          tlinlog=0
        case('lint')     ! linear time for trajectories
          tlinlog=1
        case('zero')     ! plot origin
          call set_zero
        case('outl')     ! outline of plasma grid
          call outlineon(.true.)
        case('stdl')     ! "standard" label
          call customise_label
     .     ('${HOST} ${PWD} `date` T=${TIME}',stmp)
          call clear_extralabel(stmp)
        case('mvsk')     ! movie skip
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          movie_skip=istack(iistack)
          call pop_stack(stack_req)
        case('ave')      ! number of things to average over
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          average=istack(iistack)
          call pop_stack(stack_req)
        case('t_#p')     ! number of trajectories
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call t_setnplot(istack(iistack))
          call pop_stack(stack_req)
        case('reso')     ! reset output
          call resetoutput
        case('swid')     ! part of tcl/tk interface
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call setwinid(istack(iistack))
          call pop_stack(stack_req)
        case('xyrz')     ! part of tcl/tk interface
          call getxyrz
        case('t_r+')     ! ratio true
          ratio=.true.
        case('t_r-')     ! ratio false
          ratio=.false.
        case('mcol')     ! colour index of mesh
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          mcol=istack(iistack)
          call pop_stack(stack_req)
        case('scol')     ! colour index of separatrix
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          scol=istack(iistack)
          call pop_stack(stack_req)
        case('vcol')     ! colour index of vessel
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          vcol=istack(iistack)
          call pop_stack(stack_req)
        case('ocol')     ! colour index of grid outline
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          ocol=istack(iistack)
          call pop_stack(stack_req)
        case('zcol')     ! colour index of zero lines
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_axes_color(istack(iistack))
          zcol=istack(iistack)
          call pop_stack(stack_req)
        case('hcol')     ! colour index of text header
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_header_color(istack(iistack))
          hcol=istack(iistack)
          call pop_stack(stack_req)
        case('gcol')     ! colour index of graph title
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_title_color(istack(iistack))
          gcol=istack(iistack)
          call pop_stack(stack_req)
        case('bcol')     ! colour index of bar label title
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_barlabel_color(istack(iistack))
          bcol=istack(iistack)
          call pop_stack(stack_req)
        case('xcol')     ! colour index of extralabel
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_extralabel_color(istack(iistack))
          xcol=istack(iistack)
          call pop_stack(stack_req)
        case('dcol')     ! colour index for matrix data line plots
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          dcol=istack(iistack)
          call pop_stack(stack_req)
        case('wcol')     ! colour index for wall data line plots
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          wcol=istack(iistack)
          call pop_stack(stack_req)
        case('tcol')     ! colour index for target data line plots
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          tcol=istack(iistack)
          call pop_stack(stack_req)
        case('rigz')     ! toggle set\_resignore\_regions\_to\_zero
          set_resignore_regions_to_zero=
     1     .not.set_resignore_regions_to_zero
        case('noax')     ! no axes
          call set_noax
        case('part')     ! partition colour map
          stack_req=(/0,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_partition_colors(rpop(rstack,irstack))
        case('lbsa')     ! toggle label lines
          call toggle_label_lines
        case('ytra')     ! transform the "y" coordinate for comparison to experiment
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          y_transform_l=istack(iistack)
          call pop_stack(stack_req)
        case('ecol')     ! colour index of experimental data
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          ecol=istack(iistack)
          call pop_stack(stack_req)
        case('fcol')     ! colour index for graph frame
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call set_frame_color(istack(iistack))
          fcol=istack(iistack)
          call pop_stack(stack_req)
        case('nspe','nspec')    ! new spectrometer
          stack_req=(/3,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          is1=istack(iistack-2)
          if(0.gt.is1.or.is1.ge.ns) then
            write(*,*) 'Faulty parameter is = ',is1,' for nspec!'
            return
          endif
c     Parameters on istack: <species number> <wavelength in A>
c     <delta wavelength in A> nspe
          call CalculateEmissivity(istack(iistack-2),
     &     istack(iistack-1),
     &     istack(iistack))
          call pop_stack(stack_req)
        case('spec')     ! old spectrometer
          call spect
        case('coro')     ! print coronal rates
          stack_req=(/1,1,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call callgetsa(ns,rstack(irstack),istack(iistack))
          call pop_stack(stack_req)
        case('cora')     ! coronal radiation rate
          call setcorora(ismain)
        case('coni')     ! coronal densities
          call setcoroni(ismain)
        case('shco')     ! part of tcl/tk
          continue
        case('user')     ! call a user function "b2paux"
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          call b2paux(istack(iistack))
          call pop_stack(stack_req)
        case('writ','write')    ! set writing to file
          if(.not.writing_to_file) then
            writing_to_file=.true.
            open(nwtf,file='b2plot.write')
          else
            writing_to_file=.false.
            close(nwtf)
          end if
        case('omp')      ! set jxa (OMP)
          stack_req=(/1,0,0,0,0,0,0/)
          if(bad_args(stack_req)) return
          if(-1.le.istack(iistack).and.istack(iistack).le.nx) then
            jxa = ipop(istack,iistack)
            write(*,*) 'jxa changed to ', jxa
          else
            write(*,*)
     .       'Faulty ix = ',istack(iistack),' chosen for jxa!'
            call pop_stack(stack_req)
          endif
        case('imp')      ! set jxi (IMP)
          stack_req=(/1,0,0,0,0,0,0/)
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
          stack_req=(/1,0,0,0,0,0,0/)
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

      subroutine initialize_neoclassical(neo_ic)
      integer neo_ic
      if(first_neoclassical) then
c calculate neoclassical diagnostics
        call alloc_b2mod_neoclassical(nx,ny,ns,mmx)
        first_neoclassical = .false.
      endif
      if(rerun_neoclassical) then
        call neo_diagnostics(nx, ny, ns, neo_ic, dn_neo, chi_neo,
     &   vn_neo, vw_neo)
        call neo_validity(nx, ny, ns, jsep, jxa, safetyfactor,
     &   nuab, nua, nustar, rhol, lnn, ltt, deltal, wb, deltaw,
     &   deltawn, deltawt, omega, deltab, deltapar, vtha)
        rerun_neoclassical = .false.
      endif
      return
      end subroutine initialize_neoclassical

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
      return
      end subroutine write_crazy
