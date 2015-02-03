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
*-----------------------------------------------------------------------
*.specification

      subroutine b2ytdr (ninp, nout, nx, ny, ns, nx1, ny1)
      use b2mod_types
      use b2mod_geo
      use b2mod_geo2
      use b2mod_diag
      use b2mod_wall
      use b2mod_wall2
      use b2mod_plasma
      use b2mod_sources
      use b2mod_version
      use b2mod_sputter
      use b2mod_indirect
      use b2mod_elements
      use b2mod_geo_corner
      use b2mod_input_profile
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_numerics_namelist
      use b2mod_transport_namelist
      use b2mod_transport_models
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_b2cmrc
      use b2mod_b2cmfs
      use b2mod_ppout
#ifdef USE_EIRENE
      use eirmod_braeir
#endif
      implicit none
#     include "DIMENSIONS.F"
#ifndef USE_EIRENE
#     include "COUPLE/BRAEIR.F"
#endif
*   ..input arguments (unchanged on exit)
      integer ninp(0:5), nout(0:3), nx, ny, ns, nx1, ny1
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YTDR is a post-processor for b2. It transfers the plasma state
*     to a mesh of a different size.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine is a post-processor for the b2 code. It reads the
*     result of a b2 calculation, transfers the plasma state to a mesh
*     of a different size and/or a different species bundling scheme, 
*     and writes that state to an un*formatted file. This routine does 
*     not employ any metric information other than the mesh dimensions
*     and connectivity. This routine also attempts to re-write the
*     input files to match the new mesh size and species bundling.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:3) integer array, input.
*     ninp specifies the input unit numbers as follows.
*     ninp(0): formatted; specifies the output mesh size.
*     ninp(1): un*formatted; provides the plasma state.
*     ninp(2): formatted; provides the run switches.
*     ninp(3): un*formatted; provides the (new) geometry.
*     ninp(4): un*formatted; provides the (new) atomic rates.
*     ninp(5): formatted; provides the transport switches.
*
*     nout - (0:2) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides printed output.
*     nout(1): graphical; provides the transferred plasma state.
*     nout(2): formatted; provides modified run switches.
*     nout(3): formatted; provides modified transport parameters.
*
*     nx, ny - integer, input.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*     ns - integer, input.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*     nx1, ny1 - integer.
*     nx1 and ny1 specify the mesh dimensions for the output state.
*     It will hold that 0.le.nx1 and 0.le.ny1.
*
*     (parameters in common:)
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer idum(0:12), ii(1), i1, i2, ik, inc, nnc
      integer nx_sp, ny_sp, ns_sp, nbc1, ndat, nds(ns-1)
      integer, parameter :: nkind_source=4, nkind_data=2, nkind_coeff=9
      integer boundary_namelist, neutrals_namelist, use_eirene
      integer numerics_namelist, sources_inputfile
      integer transport_namelist, transport_inputfile
      integer atomic_physics_rescale_flag
      integer ismain, ismain0, isfeedback
      integer sput_dst, sput_dst2, sput_dst3
      integer sput_phys_col, sput_phys_model, refl_model
      integer iy_nocoreExB, bc_type16_ref, bc_type20_ref, ixref, iyref
      integer ixref1, iyref1
      integer ix1, ix2, iy1, iy2
      integer ix11, ix21, iy11, iy21
      integer elm_ix1, elm_ix2
      integer mw_ix1, mw_ix2, mw_jxa, mw_jsep
      integer mw_jxa1, mw_jsep1
      integer jsep, jxa, jxi, jsep1, jxi1, jxa1, iswitch
      integer nncut1, nnreg1(0:2), periodic_bc1
      integer plate_model, plate_option
      integer scaling_ix_begin, scaling_ix_end
      integer elm_ix_begin, elm_ix_end
      integer ndim, ndm2, nfitanf, nfitend, ifail
      character lblmn*120, lblgm*120, lblrc*120, lblcp*120, ligne*120
      character*10 :: b2fgmtry_oldversion, b2fgmtry_newversion
      data b2fgmtry_oldversion / '00.000.000' /
      character tab
      real (kind=R8) ::
     1 fhe1(-1:nx1,-1:ny1,0:1), 
     2 fch1(-1:nx1,-1:ny1,0:1),
     3 fhi1(-1:nx1,-1:ny1,0:1), po1(-1:nx1,-1:ny1),
     5 ti1(-1:nx1,-1:ny1), te1(-1:nx1,-1:ny1), ne1(-1:nx1,-1:ny1), 
     6 fch1_32(-1:nx1,-1:ny1,0:1), fch1_52(-1:nx1,-1:ny1,0:1),
     7 fch1_p(-1:nx1,-1:ny1), w1d(0:ns-1)
      real (kind=R8), allocatable ::
     1 zamin1(:), zamax1(:), zn1(:), am1(:),
     2 uadia1(:,:,:,:), fna1(:,:,:,:), ua1(:,:,:), na1(:,:,:), 
     3 kinrgy1(:,:,:),
     4 sputter_yield_yt(:,:,:,:), sputter_yield2_yt(:,:,:,:), 
     . sputter_tmp(:), cftmp(:), 
     5 bundle_contr(:,:), bundle_frac(:,:),
     6 s1data(:,:,:,:), t1data(:,:,:,:), eta_hce_mult1(:),
     7 recyc1(:,:), erecyc1(:,:), rcion1(:,:), 
     8 phys_sput1(:,:), chem_sput1(:,:), 
     9 r(:), f(:), d(:), felm(:), delm(:), px(:), pf(:), pelm(:),
     * prdsa(:), pra(:), dsa(:), dsa1(:), deta(:),
     1 sr1(:), sf1(:), sd1(:), sr2(:), sf2(:), sd2(:),
     2 spr(:), spf1(:), spf2(:), s1d(:,:,:,:)
      real (kind=R8) :: ds_offset
      integer ns1, nsm, i, j, k, ix, iy, is, ib, ireg, len, nsns
      integer ib1, ib2, ifl, js, ks, kmod, itrack
      integer nspecies1, iss
      integer itrack1, itrack2, itrack3
      integer non_commensurate
      integer is0(0:100), is1(0:ns-1)
      integer is1beginrange, is1endrange, iatmbeginrange, iatmendrange
      integer iatm1(0:ns), ipls1(0:ns)
      integer, allocatable :: istmp(:), ns1data(:,:,:)
      integer, allocatable :: mlcmp1(:,:), ifl1(:)
      integer, allocatable :: b2espcr1(:), eb2spcr1(:), eb2ibcr(:)
      integer, allocatable :: e_bundle_start(:), e_bundle_end(:)
      integer, allocatable :: n1data(:,:,:), ad1spec(:,:,:), isp1(:)
      integer bcpos1(nbcd), bcstart1(nbcd), bcend1(nbcd)
      integer rcpos1(nstraid), rcstart1(nstraid), rcend1(nstraid)
      integer species_start1(nstraid), species_end1(nstraid)
      integer targsp1(nstraid,ntrack), lstrascl1(nstraid,natm)
      integer nxiso3(nnisod), nxiso4(nnisod)
      integer nyiso3(nnisod), nyiso4(nnisod)
      integer nsdata(nkind_data,nkind_source,0:DEF_NSD)
      integer ndata(nkind_data,nkind_coeff,0:DEF_NSD)
      integer addspec(DEF_NSD,nkind_coeff,0:DEF_NSD)
      integer ixa, ixb, iya, iyb
      integer nstrai1
      integer IUNIN, IUNOT, INDGRD(3), NR1ST, NP2ND, NSTSI, NLIMI
      integer INDPRO(12), NMODE, NAINI, NAINB
      integer NZADD, NVOLPR, NSURPR, NVOLPL, NSP
      integer NPLINR,NPLOTR,NPLDLR,NPLINP,NPLOTP,NPLDLP,
     .        NPLINT,NPLOTT,NPLDLT
      integer IRPTA1, IRPTE1, IRPTA2, IRPTE2
      integer NATMI, NMOLI, NIONI, NPLSI, NPLSI1, NREAC, IPLTA
      integer NATMI1
      integer IREAD, ISTS, IH, L, M
      integer IATM, IMOL, IION, IPLS, ICRM, ISTRA, NSTRAE, NSTRA1
      integer MAXRAD, NSIGVI, NSIGSI, NSIGCI, NTARGI, NPPLG
      integer NSIGV1, NSIGS1, NSIGC1
      integer NADVI, NCLVI, NALVI, NADSI, NALSI, NAOTB
      integer ISP, IPROJ, NPLI1, NPLO1, NPLI2, NPLO2
      integer NPLIN2, NPLOT2, NPLDL2
      integer NCHORI, NSNVI
      integer NDT, NINCT, NIXY, NTIN, NTEN, NIFLG, NPTC
      integer iside, numpairs
#ifndef NEW_SIZE
      integer nchsppr, nchsppri, nsrfcls
      integer NLMAX
      integer lntspec, lntsput, nlstadds, ngwties, mgwties
      character*2 txtmtmp, txtmtra
#else
      integer NTL
      integer IPOS0, IPOS1, IPOS2
      integer NPRCSF, IPRCSF
      integer MAXPOL, MAXTOR, MAXADD
      integer NSPZI, NSPZE, NSPZ1
      integer NSIGI_BGK, NSIGI_COP
#endif
      integer, allocatable :: NPOINT(:,:), ILIIN(:)
      integer, allocatable :: INDSRC(:)
      integer, allocatable :: NMASSA(:), NCHARA(:), NPRTA(:), 
     &                        ISRFA(:), ISRTA(:), ISRSA(:), NRCA(:),
     &                        IREACA(:,:),IBULKA(:,:),
     .                        ISCD1A(:,:),ISCD2A(:,:),
     .                        ISCDEA(:,:)
      integer, allocatable :: NMASSM(:), NCHARM(:), NPRTM(:), 
     &                        ISRFM(:), ISRTM(:), ISRSM(:), NRCM(:),
     &                        IREACM(:,:),IBULKM(:,:),
     .                        ISCD1M(:,:),ISCD2M(:,:),
     .                        ISCDEM(:,:)
      integer, allocatable :: NMASSI(:), NCHARI(:), NPRTI(:), NCHRGI(:),
     &                        ISRFI(:), ISRTI(:), ISRSI(:), NRCI(:),
     &                        IREACI(:,:),IBULKI(:,:),
     .                        ISCD1I(:,:),ISCD2I(:,:),
     .                        ISCDEI(:,:)
      integer, allocatable :: NMASSP(:), NCHARP(:), NPRTP(:), NCHRGP(:), 
     &                        ISRFP(:), ISRTP(:), ISRSP(:), NRCP(:),
     &                        IREACP(:,:),IBULKP(:,:),
     .                        ISCD1P(:,:),ISCD2P(:,:),
     .                        ISCDEP(:,:)
      integer, allocatable :: NMASSP1(:), NCHARP1(:), 
     &                        NPRTP1(:), NCHRGP1(:), 
     &                        ISRFP1(:), ISRTP1(:), ISRSP1(:), NRCP1(:)
      integer, allocatable :: IPLAA(:), IPLEA(:)
      integer, allocatable :: IPLAS(:), IPLES(:)
      integer, allocatable :: IFLB(:), IFLB1(:)
      integer, allocatable :: NTGPRT(:)
      integer, allocatable :: NPTS(:), NINITL(:), NEMODS(:), NAMODS(:)
      integer, allocatable :: NSPEZ(:)
      integer, allocatable :: NSRFSI(:)
      integer, allocatable :: INDIM(:,:), INSOR(:,:),
     .                        INGRDA(:,:,:), INGRDE(:,:,:)
      integer, allocatable :: NRSOR(:,:), NPSOR(:,:), NTSOR(:,:),
     .                        NBSOR(:,:), NASOR(:,:), NISOR(:,:),
     .                        ISTOR(:,:)
      integer, allocatable :: IGH(:), IIH(:), IGHW(:), IIHW(:),
     .                        IGHC(:,:), IIHC(:,:)
#ifdef NEW_SIZE
      integer, allocatable :: NFOLA(:), NGENA(:)
      integer, allocatable :: NFOLM(:), NGENM(:)
      integer, allocatable :: NFOLI(:), NGENI(:)
      integer, allocatable :: IESTMA(:,:), IBGKA(:,:)
      integer, allocatable :: IESTMM(:,:), IBGKM(:,:)
      integer, allocatable :: IESTMI(:,:), IBGKI(:,:)
      integer, allocatable :: IVLSF(:),
     .                        ISCLS(:), ISCLT(:), ISCL1(:),
     .                        ISCL2(:), ISCL3(:), ISCLB(:),
     .                        ISCLA(:)
#else
      integer, allocatable :: NSSPL(:)
      integer, allocatable :: lsrfcls(:)
      integer, allocatable :: listadds(:), listaddg(:), lgwties(:)
      integer, allocatable :: lkindp(:), lkindp1(:)
#endif
      real (kind=R8) :: rdum(0:17), RMN, RMX, ALLOC
      real (kind=R8), allocatable :: XPOL(:,:), YPOL(:,:)
      real (kind=R8), allocatable :: VOLCOR(:), VOLADD(:)
      real (kind=R8), allocatable :: RLB(:)
      real (kind=R8), allocatable :: EELECA(:,:), EBULKA(:,:),
     .           ESCD1A(:,:), ESCD2A(:,:), FREACA(:,:)
      real (kind=R8), allocatable :: EELECM(:,:), EBULKM(:,:),
     .           ESCD1M(:,:), ESCD2M(:,:), FREACM(:,:)
      real (kind=R8), allocatable :: EELECI(:,:), EBULKI(:,:),
     .           ESCD1I(:,:), ESCD2I(:,:), FREACI(:,:)
      real (kind=R8), allocatable :: EELECP(:,:), EBULKP(:,:),
     .           ESCD1P(:,:), ESCD2P(:,:),FREACP(:,:)
      real (kind=R8), allocatable :: 
     .  TIO(:), TIl(:), TI2(:), TI3(:), TI4(:), TI5(:)
      real (kind=R8), allocatable :: 
     .  DI0(:), DI1(:), DI2(:), DI3(:), DI4(:), DI5(:)
      real (kind=R8), allocatable :: 
     .  VX0(:), VX1(:), VX2(:), VX3(:), VX4(:), VX5(:)
      real (kind=R8), allocatable :: 
     .  VY0(:), VY1(:), VY2(:), VY3(:), VY4(:), VY5(:)
      real (kind=R8), allocatable :: 
     .  VZ0(:), VZ1(:), VZ2(:), VZ3(:), VZ4(:), VZ5(:)
      real (kind=R8), allocatable :: 
     .  P0(:), P1(:), P2(:), P3(:), P4(:), P5(:)
      real (kind=R8), allocatable :: DATM(:), DATM1(:)
      real (kind=R8), allocatable :: DPLS(:), DPLS1(:)
      real (kind=R8), allocatable :: FCTE(:), FCTE1(:)
      real (kind=R8), allocatable :: BMASS(:), BMASS1(:)
      real (kind=R8), allocatable :: FLUX(:)
      real (kind=R8), allocatable :: SORWGT(:,:), SORLIM(:,:),
     .                               SORIND(:,:), SOREXP(:,:),
     .                               SORIFL(:,:)
      real (kind=R8), allocatable :: SORAD1(:,:), SORAD2(:,:),
     .                               SORAD3(:,:), SORAD4(:,:),
     .                               SORAD5(:,:), SORAD6(:,:)
      real (kind=R8), allocatable :: SORENE(:), SORENI(:),
     .                               SORVDX(:), SORVDY(:), SORVDZ(:)
      real (kind=R8), allocatable :: SORCOS(:), SORMAX(:),
     .                               SORCTX(:), SORCTY(:), SORCTZ(:)
#ifdef NEW_SIZE
      real (kind=R8), allocatable :: SCALV(:)
      real (kind=R8), allocatable :: XTRIAN(:), YTRIAN(:)
#endif
      real (kind=R8) ::
     *  neutral_rescale, plate_temp, plate_thick, BoRiS, divheat,
     .  scaling_width, scaling_strength,
     .  elm_time_begin, elm_time_end, elm_time_period,
     1  zamin0(0:100), zamax0(0:100), zn0(0:100), am0(0:100),
     2  sdata(2,DEF_NXD,nkind_source,0:DEF_NSD), 
     3  tdata(3,DEF_NXD,nkind_coeff,0:DEF_NSD),
     4  rescale_ra(0:nsdecl-1), rescale_sa(0:nsdecl-1),
     5  rescale_qa(0:nsdecl-1), rescale_cx(0:nsdecl-1),
     6  rescale_rd(0:nsdecl-1), rescale_br(0:nsdecl-1)
      real (kind=R8) :: recyceir1(nstraid), volrecstart1(nstraid), 
     &  userfluxparm1(nstraid,2)
      character cnamip*80, cvalip*80, geometry_comment*80
      character*1 crcstra1(nstraid)
      character*32 hlp_frm, hlp_frp, hlp_frb
      character*80 line
      character*24 TXTSPC, TXTUNT, TEXTLA
      character*72 ZEILE, TXTMOL, TXTION, TXTPLS, ZEILEN(100)
      character*72 TITLE5, CHRTAL
      character FILNAM*8,H123*4,REAC*9,CRC*3
      character*8, allocatable :: TEXTA(:), TEXTM(:), TEXTI(:), TEXTP(:)
      character*8, allocatable :: TEXTP1(:)
      character*72, allocatable :: TXTSOU(:)
      character*10, allocatable :: species_name(:)
      logical ldum(0:60), found, ex
      logical chemsp1(nstraid), savestratum(0:nstraid)
      logical, allocatable :: is_neutral1(:), done(:)
      logical, allocatable :: savesigv(:), savesigs(:), savesigc(:)
      logical, allocatable :: NLAVRP(:), NLAVRT(:), NLSYMP(:), NLSYMT(:)
      logical, allocatable :: NLATM(:), NLMOL(:), NLION(:), NLPLS(:)
      logical, allocatable :: NLPNT(:), NLLNE(:), NLSRF(:),
     .                        NLVOL(:), NLCNS(:),
     .                        NLPRCA(:), NLPRCM(:), NLPRCI(:), 
     .                        NLPRCA1(:)
      logical, allocatable :: TRCSRC(:), TRCSR1(:)
#ifndef NEW_SIZE
      logical, allocatable :: nlpbls(:), nlpbl1(:)
#endif
      logical todo, same_bundling, in_bundle
      logical sputter_namelist_used, wall_namelist_used
      logical region_flags(NNREGMAX,nkind_coeff), no_pflux
      logical poloidal_scaling, prdone, clash(ns-1), shared(0:ns-1)
      logical track_chem_sput1(ntrack)
*   ..procedures
      integer lnblnk, convert_position, len_of_digits
      integer conformal_convert, convert_eirene, b2idez
      real (kind=R8) :: cr, cz, cr1, cz1
      logical open_file
      external subini, subend, xertst, cfruch, cfrure, cfwuin,
     &  cfwuch, cfwure, ipgeti, lnblnk, convert_position,
     &  len_of_digits, open_file, get_jsep, conformal_convert,
     &  convert_eirene, b2idez
      external b2rups, b2rurc, b2wups, b2xtps, b2rflb, b2ytgm
      external b2ahx0, b2rfcp, b2wfcp, b2xvcp, b2rugm, cfverr
      external E01BEF, E01BFF
      intrinsic nint, min0, max, index, iabs, sqrt
*   ..initialisation
      data non_commensurate/0/
      data boundary_namelist/0/, neutrals_namelist/0/, use_eirene/0/
      data numerics_namelist/0/, sources_inputfile/0/
      data transport_namelist/0/, transport_inputfile/0/
      data atomic_physics_rescale_flag/0/
      data sput_dst/-1/, sput_dst2/-1/, sput_dst3/-1/
      data itrack1/-1/, itrack2/-1/, itrack3/-1/
      data ismain/1/, ismain0/0/, isfeedback/0/
      data sput_phys_col/3/, sput_phys_model/1/, refl_model/1/
      data iy_nocoreExB/-2/, bc_type16_ref/1/, bc_type20_ref/1/
      data ixref/-1/, iyref/-1/, mw_jxa/-1/, mw_jsep/-1/
      data jxa/-2/, jxi/-2/
      data ix1/-1/, ix2/-1/, iy1/-1/, iy2/-1/
      data elm_ix1/-1/, elm_ix2/-1/, mw_ix1/-1/, mw_ix2/-1/
      data plate_temp/300.0_R8/, plate_thick/0.00_R8/
      data plate_model/0/, plate_option/3/
      data BoRiS/0.0_R8/
*  ..namelists
      namelist /sputter_save/ nx_sp, ny_sp, ns_sp,
     *   sputter_yield, sputter_yield2
      namelist /profile/ nsdata, sdata, divheat
      namelist /transport/ ndata, tdata, addspec,
     1   transport_filename, transport_time_mod, transport_time_switch,
     2   scaling_ix_begin, scaling_ix_end, scaling_strength,
     3   region_flags, no_pflux, poloidal_scaling, scaling_width,
     4   elm_time_begin, elm_time_end, elm_time_period,
     5   elm_ix_begin, elm_ix_end
      namelist /atomic_physics_rescale/
     1 rescale_sa, rescale_ra, rescale_qa, rescale_cx,
     3 rescale_rd, rescale_br

      cr(ix,iy)=0.25_R8*
     1 (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
      cz(ix,iy)=0.25_R8*
     1 (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))
      cr1(ix,iy)=0.25_R8*
     1 (crx1(ix,iy,0)+crx1(ix,iy,1)+crx1(ix,iy,2)+crx1(ix,iy,3))
      cz1(ix,iy)=0.25_R8*
     1 (cry1(ix,iy,0)+cry1(ix,iy,1)+cry1(ix,iy,2)+cry1(ix,iy,3))

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The workspace array w is used as a stack for storage of dynamic
*     real*8 arrays. An array 'name' in w(0:*) will be stored starting at
*     location w(l'name').
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     /b2cmpa/ specifies basic physics parameters.
*      me, mp, ev, qe, eps0, mu0 - real*8 constants.
*     These parameters are included with /b2cmpa/.
*       me = electron mass (kg)
*       mp = proton mass (kg)
*       ev = 1 electron-Volt (Joule)
*       qe = elementary charge (Coulomb)
*       eps0 = vacuum dielectric permittivity (SI)
*       mu0 = vacuum magnetic permeability (SI)
*      nsdecl - integer constant.
*     nsdecl is an upper bound for the number of species (ns), and is
*     used in dimension specifications.
*      zamin, zamax, zn, am - (0:nsdecl-1) real*8 array.
*     For (is) in (0:ns-1), zn(is) specifies the nuclear charge of
*     species (is), am(is) specifies its atomic mass, and zamin(is) and
*     zamax(is) specify its atomic charge range. All are dimensionless.
*     It will hold that 0.le.zamin(is).le.zamax(is), 0.lt.zn(is) and
*     0.lt.am(is). It is required that neutral species not be bundled
*     with charged species.
*     For ns.le.is.lt.nsdecl, zamin(is), zamax(is), zn(is), am(is) 
*     are undefined.
*
*     Description of some local variables (including arrays in w(0:*):
*
*     na - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), na(ix,iy,is) specifies
*     the density of atomic species (is) on the (ix,iy) cell.
*     It will hold that 0.lt.na(,,).
*
*     ne - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ne(ix,iy) specifies
*     the electron density on the (ix,iy) cell.
*     It will hold that 0.lt.ne(,).
*
*     ua - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), ua(ix,iy,is) specifies
*     the parallel velocity of atomic species (is) on the (ix,iy) cell.
*
*     uadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), uadia(ix,iy,0,is) specifies
*     the diamagnetic drift velocity of species (is) on the (ix,iy) cell, and
*     uadia(ix,iy,1,is) the radial drift velocity.
* 
*     te - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), te(ix,iy) specifies the electron
*     temperature on the (ix,iy) cell. It will hold that 0.lt.te(,).
*
*     ti - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ti(ix,iy) specifies the all atom
*     temperature on the (ix,iy) cell. It will hold that 0.lt.ti(,).
*
*     po - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), po(ix,iy) specifies the electric
*     potential on the (ix,iy) cell.
*
*     fna - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let fnax(,,)=fna(,,0,) and fnay(,,)=fna(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), fnax(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fnax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fnay(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fnay(-1:nx,iy,0:ns-1) will hold 0.
*
*     fhe - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhex(,)=fhe(,,0) and fhey(,)=fhe(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhex(ix,iy) specifies the electron
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhex(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhey(ix,iy) specifies the electron
*     heat flux through the face between the (ix,iy-1) and (ix,iy)
*     cells. For iy.eq.-1, fhey(-1:nx,iy) will hold 0.
*
*     fhi - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhix(,)=fhi(,,0) and fhiy(,)=fhi(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhix(ix,iy) specifies the all atom
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhiy(ix,iy) specifies the all atom
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For iy.eq.-1, fhiy(-1:nx,iy) will hold 0.
*
*     fch - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchx(,)=fch(,,0) and fchy(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchx(ix,iy) specifies the electric
*     current through the face between the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchy(ix,iy) specifies the electric
*     current through the face between the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fchy(-1:nx,iy) will hold 0.
*
*     na1, ne1, ua1, uadia1, te1, ti1, po1, fna1, fhe1, fhi1, fch1 - real*8 arrays.
*     These arrays hold the state on the (-1:nx1,-1:ny1) mesh.
*
*     lblmn - character*120
*     lblmn identifies the instance of the program employed for the main
*     calculation.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ytdr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.
     &             1.le.ninp(3), 'faulty argument ninp')
      call xertst (1.le.nout(0).and.1.le.nout(1).and.1.le.nout(2),
     &  'faulty argument nout')
*   ..test dimensions
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')
      call xertst (0.le.nx1.and.0.le.ny1, 'faulty argument nx1, ny1')
      write(*,*) 'nx, nx1 = ', nx, nx1
      write(*,*) 'ny, ny1 = ', ny, ny1
*   ..read header for geometry
      call cfruch (ninp(3), 120, lblgm, 'label')
      call cfruin (ninp(3), 1, idum, 'isymm')
      isymm = idum(0)
*     (isymm is now stored in the b2cmfs common)               ! xpb 10.08.00
*   ..read new geometry and metric
      call b2ytgm (ninp(3), nx1, ny1, nncut, nnreg)
      nncut1 = nncut
      nnreg1 = nnreg
      periodic_bc1 = periodic_bc
      call b2yt_jsep (nx1, ny1, jxi1, jxa1, jsep1)
*   ..read and echo code internal parameters from b2yt.dat
      write (nout(0),'(/2x,a)') 
     . 'non-default internal parameters from b2yt.dat:'
    3 continue
       read (ninp(0),*,end=4,err=91) cnamip, cvalip
       write (nout(0),'(4x,a,2x,a)') cnamip, cvalip
       call ipsetc (cnamip, cvalip)
       goto 3
    4  continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'

* ..computation
*   ..read plasma state file label
      call cfruch (ninp(1), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(1), b2fstate_version, ns, zamin, zamax, zn, am)
*   ..read plasma state
      call b2rups (ninp(1), nx, ny, ns, ne, na, ua, uadia, te,
     &  ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     1  b2fstate_version)
*   ..write label
      write (nout(0),'(/2x,a)') 'label for this calculation:'
      write (nout(0),'(/4x,a/4x,a)') lblmn(1:60), lblmn(61:120)
*   ..write mesh dimensions
      write (nout(0),101) nx1, ny1
*   ..check if new bundling or species are required and allocate ns-dependent arrays
      ns1=ns
      call ipgeti ('b2ytdr_ns', ns1)
      call xertst (0.lt.ns1, 'faulty input ns1')
      call xertst (DEF_NSD.ge.ns1, 
     . 'increase size of DEF_NSD in DIMENSIONS.F')
      call ipgeti ('b2ytdr_non_commensurate', non_commensurate)
      if (non_commensurate.eq.1) then
*   ..check for and open old geometry file
        call xertst (open_file(70,'b2fgmtry.old',''),
     .   'Non-commensurate translation requires old b2fgmtry.old file')
*   ..obtain old geometry version number
        call cfverr (70, b2fgmtry_oldversion)
        if (b2fgmtry_oldversion.ne.b2fgmtry_version) then
          b2fgmtry_newversion=b2fgmtry_version
          b2fgmtry_version=b2fgmtry_oldversion
        endif
*   ..obtain nx, ny from old geometry file
        call cfruin (70, 2, idum, 'nx,ny')
        call xertst (nx.eq.idum(0), 
     &   'inconsitent nx definition between b2fgmtry.old and b2yt.dat!')
        call xertst (ny.eq.idum(1),
     &   'inconsitent ny definition between b2fgmtry.old and b2yt.dat!')
*   ..read header for old geometry
        call cfruch (70, 120, lblgm, 'label')
        call cfruin (70, 1, idum, 'isymm')
*   ..read old geometry and metric
        call b2rugm (70, nx, ny, crx, cry, fpsi,
     &   ffbz, bb, vol, hx, hy, qz, qc, gs, pbs, wbbl, wbbr, wbbv)
        call xertst (nnreg(0).eq.nnreg1(0) .and. nncut.eq.nncut1 .and.
     &   periodic_bc.eq.periodic_bc1,
     &   'Non-commensurate translation requires the same topology!')
        if (nncut.eq.2) call xertst(
     .   (topcut(1).eq.topcut(2) .and. topcut1(1).eq.topcut1(2)) .or.
     .   (topcut(1).lt.topcut(2) .and. topcut1(1).lt.topcut1(2)) .or.
     .   (topcut(1).gt.topcut(2) .and. topcut1(1).gt.topcut1(2)),
     .   'Non-commensurate translation requires the same DN topology!')
        call get_jsep(nx,ny,jxi,jxa,jsep)
        close(70)
        if (b2fgmtry_oldversion.ne.b2fgmtry_newversion) then
          b2fgmtry_version=b2fgmtry_oldversion
        endif
      endif
      nsm=max(ns,ns1)
      allocate(am1(0:ns1-1))
      allocate(zn1(0:ns1-1))
      allocate(zamin1(0:ns1-1))
      allocate(zamax1(0:ns1-1))
      allocate(is_neutral1(0:ns1-1))
      allocate(na1(-1:nx1,-1:ny1,0:ns1-1))
      allocate(ua1(-1:nx1,-1:ny1,0:ns1-1))
      allocate(kinrgy1(-1:nx1,-1:ny1,0:ns1-1))
      allocate(fna1(-1:nx1,-1:ny1,0:1,0:ns1-1))
      allocate(uadia1(-1:nx1,-1:ny1,0:1,0:ns1-1))
      allocate(sputter_yield_yt(-1:nx1,-1:ny1,0:ns1-1,1:2))
      allocate(sputter_yield2_yt(-1:nx1,-1:ny1,0:ns1-1,1:2))
      allocate(bundle_frac(0:ns1-1,0:ns-1))
      allocate(bundle_contr(0:ns-1,0:100))
      allocate(b2espcr1(0:ns1-1))
      allocate(eb2spcr1(DEF_NATM))
      allocate(eb2ibcr(DEF_NPLS))
      allocate(isp1(DEF_NATM))
      allocate(ifl1(DEF_NPLS))
      allocate(e_bundle_start(0:ns1-1))
      allocate(e_bundle_end(0:ns1-1))
      allocate(sputter_tmp(0:nsm-1))
      allocate(cftmp(0:nsm-1))
      allocate(istmp(0:nsm-1))
      allocate(recyc1(0:ns1-1,DEF_NSTRA))
      allocate(erecyc1(0:ns1-1,DEF_NSTRA))
      allocate(rcion1(0:ns1-1,DEF_NSTRA))
      allocate(phys_sput1(0:ns1-1,DEF_NSTRA))
      allocate(chem_sput1(0:ns1-1,DEF_NSTRA))
      allocate(mlcmp1(DEF_NATM,DEF_NMOL))
      bundle_frac=0.0_R8
      bundle_contr=0.0_R8
      eb2spcr1=-1
      b2espcr1=-1
*   ..read new atomic rate data
      call cfruch (ninp(4), 120, lblrc, 'label')
      call alloc_b2mod_ppout(nx1,ny1,ns1)
      call read_b2mod_ppout('rates', ninp(4), lblrc)
      call b2rurc (ninp(4), b2frates_version)
*   ..save and test rtns, rtzmin, rtzmax, rtzn
      call xertst (ns1.eq.rtns, 'faulty input rtns--wrong rate table?')
      zamin1(0:ns1-1) = rtzmin(0:ns1-1)
      zamax1(0:ns1-1) = rtzmax(0:ns1-1)
      zn1(0:ns1-1) = rtzn(0:ns1-1)
      am1(0:ns1-1) = rtzn(0:ns1-1)*2.0_R8
      is_neutral1(0:ns1-1) = 
     .   nint(zamin1(0:ns1-1)).eq.0 .and. nint(zamax1(0:ns1-1)).eq.0
      same_bundling = ns1.eq.ns
      is = 0
      do while (same_bundling .and. is.lt.ns1)
        same_bundling = same_bundling .and.
     &   nint(zamin1(is)).eq.nint(zamin(is)) .and.
     &   nint(zamax1(is)).eq.nint(zamax(is)) .and.
     &   nint(zn1(is)).eq.nint(zn(is))
        is = is + 1
      end do
      nspecies = 1
      ib = 1
      do while (ib.lt.ns)
        if (nint(zn(ib)).ne.nint(zn(ib-1)) .or.
     &      nint(am(ib)).ne.nint(am(ib-1))) nspecies=nspecies+1
        ib = ib + 1
      end do
      if (same_bundling) then
        do is = 0, ns-1
          is1(is) = is
          am1(is) = am(is)
          zn1(is) = zn(is)
          zamin1(is) = zamin(is)
          zamax1(is) = zamax(is)
          bundle_frac(is,is) = 1.0_R8
        enddo
        nspecies1 = nspecies
      else 
        call de_bundle (ns, zn, am, zamin, zamax,
     *   is0, bundle_contr, zn0, am0, zamin0, zamax0)
        call re_bundle (ns, is0, bundle_contr, zn0, am0, zamin0, zamax0,
     *   ns1, is1, bundle_frac, zn1, am1, zamin1, zamax1)
        nspecies1 = 1
        ib = 1
        do while (ib.lt.ns1)
          if (nint(zn1(ib)).ne.nint(zn1(ib-1)) .or.
     &        nint(am1(ib)).ne.nint(am1(ib-1))) nspecies1=nspecies1+1
          ib = ib + 1
        end do
      endif
      ifl = 0
      do is = 0, ns1-1
        if (.not.is_neutral1(is)) then
          ifl = ifl + 1
          ifl1(is) = ifl
        end if
      end do
*   ..compute new state
      call b2xtps (nx, ny, ns, nx1, ny1, ns1, nncutmax, 
     &    non_commensurate,
     &    nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &    leftcut, leftcut1, rightcut, rightcut1,
     &    bottomcut, bottomcut1, topcut, topcut1,
     &    hx, hy, qc, hx1, hy1, qc1, 
     &    bundle_frac, bottomiy1, ne, na, ua, uadia, te, ti, 
     &    po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &    ne1, na1, ua1, uadia1, te1, ti1, po1, fna1, fhe1, fhi1,
     &    fch1, fch1_32, fch1_52, fch1_p, kinrgy1)
*   ..rescale transformed state (if necessary)
      neutral_rescale=1.0_R8
      call ipgetr ('b2ytdr_rescale_neutrals', neutral_rescale)
      call xertst (0.0_R8.lt.neutral_rescale, 
     & 'faulty input parameter neutral_rescale')
      call alloc_b2mod_elements(ns,zn,am)
      if(neutral_rescale.ne.1.0_R8) then
        write(*,*) 'Rescaling neutrals by ',neutral_rescale
        do is=0,ns1-1
          if(is_neutral1(is)) then
            na1(:,:,is)=na1(:,:,is)*neutral_rescale
            fna1(:,:,:,is)=fna1(:,:,:,is)*neutral_rescale
            kinrgy1(:,:,is)=kinrgy1(:,:,is)*neutral_rescale
          endif
        enddo
      endif

*   ..write header for transformed state
      write(*,*) 'ns, ns1 = ', ns, ns1
      idum(0) = nx1
      idum(1) = ny1
      idum(2) = ns1
      call cfverw (nout(1), newversion)
      call cfwuin (nout(1), 3, idum, 'nx,ny,ns')
      call cfwuch (nout(1), 120, lblmn, 'label')
*   ..save /b2cmpa/
      call b2wuzd (nout(1), newversion, ns1, zamin1, zamax1, zn1, am1)
*   ..save transformed state
      call b2wups (nout(1), nx1, ny1, ns1, 
     &  ne1, na1, ua1, uadia1, te1, ti1, po1, 
     &  fna1, fhe1, fhi1, fch1, fch1_32, fch1_52, fch1_p, kinrgy1)

* ..check for need to transform other input files: 
* ..b2mn.dat, b2.*.parameters, input.dat, b2ah.dat

*   ..now we treat b2mn.dat
*   ..rewrite b2mn.dat into b2mn.dat.2
*   ..read in label and revised parameters in physics common
      call b2rflb (ninp(2), lblmn(61:120))
*   ..take into account the transport parameter arrays
      write(nout(2),'(a)') 
     . '*label          (lblmn, character*60; free format)' 
      write(nout(2),*) '''',trim(lblmn(61:120)),''''
      call b2ytcp (ninp(2), nout(2), ns, ns1,
     . zamin1, zamax1, zn1, am1, is_neutral1, bundle_frac)
*   ..read and echo code internal parameters from b2mn.dat
      write (nout(0),'(/2x,a)') 
     . 'non-default internal parameters from b2mn.dat:'
    1 continue
       read (ninp(2),*,end=2,err=91) cnamip, cvalip
       write (nout(0),'(4x,a,2x,a)') cnamip, cvalip
       call ipsetc (cnamip, cvalip)
       goto 1
    2 continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'
*   ..check parameters likely to need modification
      call ipgeti ('b2stbc_boundary_namelist', boundary_namelist)
      call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)
      call ipgeti ('b2srdt_numerics_namelist', numerics_namelist)
      call ipgeti ('b2tqna_transport_namelist', transport_namelist)
      call ipgeti ('b2tqna_inputfile', transport_inputfile)
      call ipgeti ('b2sral_inputfile', sources_inputfile)
      call ipgeti ('b2mndr_atomic_physics_rescale', 
     .              atomic_physics_rescale_flag)
      call ipgeti ('b2mndr_eirene', use_eirene)
      call ipgeti ('b2mndr_ismain', ismain)
      call ipgeti ('b2mwti_ismain0', ismain0)
      call ipgeti ('b2stbr_sput_dst', sput_dst)
      call xertst (ns.gt.sput_dst, 'faulty input sput_dst')
      call ipgeti ('b2stbr_sput_dst2', sput_dst2)
      call ipgeti ('b2stbr_sput_dst3', sput_dst3)
      call xertst (sput_dst2.lt.ns, 'faulty b2stbr_sput_dst2')
      call xertst (sput_dst3.lt.ns, 'faulty b2stbr_sput_dst3')
      call ipgeti ('b2news_iy_nocoreExB', iy_nocoreExB)
      call ipgeti ('b2stbc_type16_ref', bc_type16_ref)
      call ipgeti ('b2stbc_type20_ref', bc_type20_ref)
      call ipgeti ('b2stbc_isfeedback', isfeedback)
      call ipgeti ('b2mwti_jxi', jxi)
      call ipgeti ('b2mwti_jxa', jxa)
      call ipgeti ('b2tqna_ixref', ixref)
      call ipgeti ('set_transport_ixref', ixref)
      call ipgeti ('set_transport_iyref', iyref)
      call ipgeti ('b2tqna_iy1', iy1)
      call ipgeti ('b2tqna_iy2', iy2)
      call ipgeti ('b2tqna_ix1', ix1)
      call ipgeti ('b2tqna_ix2', ix2)
      call ipgeti ('b2tqna_elm_ix1', elm_ix1)
      call ipgeti ('b2tqna_elm_ix2', elm_ix2)
      call ipgeti ('b2tqna_mw_ix1', mw_ix1)
      call ipgeti ('b2tqna_mw_ix2', mw_ix2)
      call ipgeti ('b2tqna_mw_jxa', mw_jxa)
      call ipgeti ('b2tqna_mw_jsep', mw_jsep)
*   ..we do the allocation calls now so that
*   ..the default setting correspond to the old
*   ..bundling values
      if (boundary_namelist.eq.1) 
     > call alloc_b2mod_boundary(ns, nbcd, nnisod)
      inquire(file='../b2.sputter_save.parameters',
     1 exist=sputter_namelist_used)
      inquire(file='../b2.wall_save.parameters',
     1 exist=wall_namelist_used)
      if (sputter_namelist_used.or.wall_namelist_used) 
     > call alloc_b2mod_sputter(nx,ny,ns,sput_dst,sput_dst2,sput_dst3)
      if (transport_namelist.eq.1) call alloc_b2mod_transport(ns)
      if (numerics_namelist.ge.1) call alloc_b2mod_numerics(ns)
      inquire(file='../b2.neutrals.parameters',
     1 exist=todo)
      if (neutrals_namelist.eq.1.or.use_eirene.ne.0.or.todo) then
        call alloc_b2mod_neutrals(ns, nstraid, sput_dst)
        call read_neutrals_namelist(nx, ny, ns, .false.)
        if (nstrai.gt.nstraid) stop 'Increase DEF_NSTRA in DIMENSIONS.F'
        call alloc_b2mod_diag(nx1,ny1,ns,nnreg(0),nstrai,natm)
      endif
*   ..modify b2mn.dat parameters
      ismain = is1(ismain)
      ismain0 = is1(ismain0)
      isfeedback = is1(isfeedback)
      itrack = 0
      if (sput_dst.ge.0) then
        if (is1(sput_dst).gt.0) then
          sput_dst = is1(sput_dst)
          itrack1 = 1
          itrack = 1
        else
          sput_dst = -1
        end if
      end if
      if (sput_dst2.ge.0) then
        if (is1(sput_dst2).ge.0) then
          if (sput_dst.lt.0) then
            sput_dst = is1(sput_dst2)
            itrack1 = 2
            itrack = 1
          else
            sput_dst2 = is1(sput_dst2)
            itrack2 = 2
            itrack = 2
          endif
        else
          sput_dst2 = -1
        endif
      endif
      if (sput_dst3.ge.0) then
        if (is1(sput_dst3).ge.0) then
          if (sput_dst.lt.0 .and. sput_dst2.lt.0) then
            sput_dst = is1(sput_dst3)
            itrack1 = 3
            sput_dst2 = -1
            sput_dst3 = -1
            itrack = 1
          else if (sput_dst2.lt.0) then
            sput_dst2 = is1(sput_dst3)
            itrack2 = 3
            sput_dst3 = -1
            itrack = 2
          else
            sput_dst3 = is1(sput_dst3)
            itrack3 = 3
            itrack = 3
          endif
        else
          sput_dst3 = -1
        endif
      endif
      if (iy_nocoreExB.gt.-2) iy_nocoreExB = conformal_convert
     .  (non_commensurate,2,jxa,iy_nocoreExB,nx,ny,nx1,ny1,2)
      bc_type16_ref = 2 + conformal_convert
     .  (non_commensurate,2,jxa,bc_type16_ref-2,nx,ny,nx1,ny1,2)
      bc_type20_ref = 2 + conformal_convert
     .  (non_commensurate,2,jxa,bc_type20_ref-2,nx,ny,nx1,ny1,2)
      if (ixref.ge.0) then
        ixref1=conformal_convert
     .   (non_commensurate,1,ixref,iyref,nx,ny,nx1,ny1,4)
      else
        ixref1=jxa1
      endif
      if (iyref.ge.0) then
        iyref1=conformal_convert
     .   (non_commensurate,2,ixref,iyref,nx,ny,nx1,ny1,4)
      else
        iyref1=jsep1
      endif
      iy11 = conformal_convert
     .   (non_commensurate,2,ix1,iy1,nx,ny,nx1,ny1,0)
      iy21 = conformal_convert
     .   (non_commensurate,2,ix2,iy2,nx,ny,nx1,ny1,2)
      ix11 = conformal_convert
     .   (non_commensurate,1,ix1,iy1,nx,ny,nx1,ny1,1)
      ix21 = conformal_convert
     .   (non_commensurate,1,ix2,iy2,nx,ny,nx1,ny1,3)
      elm_ix1 =conformal_convert
     .   (non_commensurate,1,elm_ix1,iy1,nx,ny,nx1,ny1,1)
      elm_ix2 = conformal_convert
     .   (non_commensurate,1,elm_ix2,iy2,nx,ny,nx1,ny1,3)
      mw_ix1 = conformal_convert
     .   (non_commensurate,1,mw_ix1,mw_jsep,nx,ny,nx1,ny1,1)
      mw_ix2 = conformal_convert
     .   (non_commensurate,1,mw_ix2,mw_jsep,nx,ny,nx1,ny1,3)
      if (mw_jxa.ge.0) then
        mw_jxa1= conformal_convert(non_commensurate,1,mw_jxa,mw_jsep,
     .   nx,ny,nx1,ny1,4)
      else
        mw_jxa1=jxa1
      endif
      if (mw_jsep.ge.0) then
        mw_jsep1= conformal_convert(non_commensurate,2,mw_jxa,mw_jsep,
     .   nx,ny,nx1,ny1,2)
      else
        mw_jsep1=jsep1
      endif
      if (jxa.ne.-2) then
        jxa1=conformal_convert(non_commensurate,1,jxa,jsep,
     .   nx,ny,nx1,ny1,4)
      endif
      if (jxi.ne.-2) then
        jxi1=conformal_convert(non_commensurate,1,jxi,jsep,
     .   nx,ny,nx1,ny1,4)
      endif
      tab = char(9)
      rewind(ninp(2))
      ligne = ' '
      do while (ligne(1:7).ne.'*endphy')
        read (ninp(2),*,end=6,err=91) ligne
      enddo
    5 continue
      read (ninp(2),*,end=6,err=91) cnamip, cvalip
      todo = .true.
      SWITCHES: select case (cnamip)
        case('b2mndr_ismain')
          iswitch = ismain
        case('b2mwti_ismain0')
          iswitch = ismain0
        case('b2stbr_sput_dst')
          iswitch = sput_dst
        case('b2stbr_sput_dst2')
          iswitch = sput_dst2
        case('b2stbr_sput_dst3')
          iswitch = sput_dst3
        case('b2stbc_isfeedback')
          iswitch = isfeedback
        case('b2news_iy_nocoreExB')
          iswitch = iy_nocoreExB
        case('b2stbc_type16_ref')
          iswitch = bc_type16_ref
        case('b2stbc_type20_ref')
          iswitch = bc_type20_ref
        case('b2tqna_ixref','set_transport_ixref')
          iswitch = ixref1
        case('set_transport_iyref')
          iswitch = iyref1
        case('b2tqna_iy1')
          iswitch = iy11
        case('b2tqna_iy2')
          iswitch = iy21
        case('b2tqna_ix1')
          iswitch = ix11
        case('b2tqna_ix2')
          iswitch = ix21
        case('b2tqna_elm_ix1')
          iswitch = elm_ix1
        case('b2tqna_elm_ix2')
          iswitch = elm_ix2
        case('b2tqna_mw_ix1')
          iswitch = mw_ix1
        case('b2tqna_mw_ix2')
          iswitch = mw_ix2
        case('b2tqna_mw_jxa')
          iswitch = mw_jxa1
        case('b2tqna_mw_jsep')
          iswitch = mw_jsep1
        case('b2mwti_jxa')
          iswitch = jxa1
        case('b2mwti_jxi')
          iswitch = jxi1
        case default SWITCHES
          todo = .false.
      end select SWITCHES
      if (todo) then
        write(hlp_frm,'(a,i1,a)') '(i',len_of_digits(iswitch),')'
        write(cvalip,hlp_frm) iswitch
      end if
      write(nout(2),'(a1,a,a,a,a1)') 
     .  '''',trim(cnamip),''''//tab//'''',trim(cvalip),''''
      goto 5
    6 continue

*  ..rewrite boundary namelist
      inquire(file='../b2.boundary.parameters',exist=todo)
      if (boundary_namelist.eq.1.or.todo) then
        call read_b2mod_boundary_namelist(nx, ny, ns, .false.)
        if(nbc.gt.nbcd) stop 'increase size in b2ytdr for nbcd'
        if(nniso.gt.nnisod) stop 'increase size in b2ytdr for nnisod'  
        nbc1 = nbc   
        do ib = 1, nbc
          call convert_bundle_i
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,bccon(0,ib))
          call convert_bundle_i
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,bcmom(0,ib))
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,conpar(0,ib,1))
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,conpar(0,ib,2))
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,conpar(0,ib,3))
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,mompar(0,ib,1))
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,mompar(0,ib,2))
          if (bcchar(ib).eq.'S') then
            bcpos1(ib) = conformal_convert(non_commensurate,
     .       2,bcstart(ib),bcpos(ib),nx,ny,nx1,ny1,2)
            bcstart1(ib) = conformal_convert(non_commensurate,
     .       1,bcstart(ib),bcpos(ib),nx,ny,nx1,ny1,1)
            bcend1(ib) = conformal_convert(non_commensurate,
     .       1,bcend(ib),bcpos(ib),nx,ny,nx1,ny1,3)
          else if (bcchar(ib).eq.'N') then
            bcpos1(ib) = conformal_convert(non_commensurate,
     .       2,bcstart(ib),bcpos(ib),nx,ny,nx1,ny1,0)
            bcstart1(ib) = conformal_convert(non_commensurate,
     .       1,bcstart(ib),bcpos(ib),nx,ny,nx1,ny1,1)
            bcend1(ib) = conformal_convert(non_commensurate,
     .       1,bcend(ib),bcpos(ib),nx,ny,nx1,ny1,3)
          else if (bcchar(ib).eq.'X') then
            bcpos1(ib) = conformal_convert(non_commensurate,
     .       2,bcstart(ib),bcpos(ib),nx,ny,nx1,ny1,4)
            bcstart1(ib) = conformal_convert(non_commensurate,
     .       1,bcstart(ib),bcpos(ib),nx,ny,nx1,ny1,1)
            bcend1(ib) = conformal_convert(non_commensurate,
     .       1,bcend(ib),bcpos(ib),nx,ny,nx1,ny1,3)
            ib1 = conformal_convert(non_commensurate,
     .       2,bcstart(ib),bcpos(ib),nx,ny,nx1,ny1,0)
            ib2 = conformal_convert(non_commensurate,
     .       2,bcstart(ib),bcend(ib),nx,ny,nx1,ny1,2)
            if (ib2.gt.ib1) then  ! must duplicate the X b.c.
              nbc1 = nbc1 + ib2-ib1
              if(nbc1.gt.nbcd) stop 'increase size in b2ytdr for nbcd'
              do is = 1, ib2-ib1
                bcchar(is) = bcchar(ib)
                bcpos1(is) = bcpos(ib)+is
                bcstart1(is) = bcstart(ib)
                bcend1(is) = bcend(ib)
                bcene(is) = bcene(ib)
                bceni(is) = bceni(ib)
                bcpot(is) = bcpot(ib)
                bccon(:,is) = bccon(:,ib)
                bcmom(:,is) = bcmom(:,ib)
                enepar(is,:) = enepar(ib,:)
                enipar(is,:) = enipar(ib,:)
                potpar(is,:) = potpar(ib,:)
                conpar(:,is,:) = conpar(:,ib,:)
                mompar(:,is,:) = mompar(:,ib,:)
                lcbs(is) = lcbs(ib)
              end do
            end if
          else if (bcchar(ib).eq.'W') then
            bcpos1(ib) = conformal_convert(non_commensurate,
     .       1,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,3)
            bcstart1(ib) = conformal_convert(non_commensurate,
     .       2,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,0)
            bcend1(ib) = conformal_convert(non_commensurate,
     .       2,bcpos(ib),bcend(ib),nx,ny,nx1,ny1,2)
          else if (bcchar(ib).eq.'E') then
            bcpos1(ib) = conformal_convert(non_commensurate,
     .       1,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,1)
            bcstart1(ib) = conformal_convert(non_commensurate,
     .       2,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,0)
            bcend1(ib) = conformal_convert(non_commensurate,
     .       2,bcpos(ib),bcend(ib),nx,ny,nx1,ny1,2)
          else if (bcchar(ib).eq.'Y') then
            bcpos1(ib) = conformal_convert(non_commensurate,
     .       1,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,4)
            bcstart1(ib) = conformal_convert(non_commensurate,
     .       2,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,0)
            bcend1(ib) = conformal_convert(non_commensurate,
     .       2,bcpos(ib),bcend(ib),nx,ny,nx1,ny1,2)
            ib1 = conformal_convert(non_commensurate,
     .       1,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,1)
            ib2 = conformal_convert(non_commensurate,
     .       1,bcpos(ib),bcend(ib),nx,ny,nx1,ny1,3)
            if (ib2.gt.ib1) then  ! must duplicate the Y b.c.
              nbc1 = nbc1 + ib2-ib1
              if(nbc1.gt.nbcd) stop 'increase size in b2ytdr for nbcd'
              do is = 1, ib2-ib1
                bcchar(is) = bcchar(ib)
                bcpos1(is) = bcpos(ib)+is
                bcstart1(is) = bcstart(ib)
                bcend1(is) = bcend(ib)
                bcene(is) = bcene(ib)
                bceni(is) = bceni(ib)
                bcpot(is) = bcpot(ib)
                bccon(:,is) = bccon(:,ib)
                bcmom(:,is) = bcmom(:,ib)
                enepar(is,:) = enepar(ib,:)
                enipar(is,:) = enipar(ib,:)
                potpar(is,:) = potpar(ib,:)
                conpar(:,is,:) = conpar(:,ib,:)
                mompar(:,is,:) = mompar(:,ib,:)
                lcbs(is) = lcbs(ib)
              end do
            end if
          endif
        enddo
        do ib = 1, nniso
          nxiso3(ib) = conformal_convert(non_commensurate,
     .       1,nxiso1(ib),nyiso1(ib),nx,ny,nx1,ny1,1)
          nxiso4(ib) = conformal_convert(non_commensurate,
     .       1,nxiso2(ib),nyiso2(ib),nx,ny,nx1,ny1,3)
          nyiso3(ib) = conformal_convert(non_commensurate,
     .       2,nxiso1(ib),nyiso1(ib),nx,ny,nx1,ny1,0)
          nyiso4(ib) = conformal_convert(non_commensurate,
     .       2,nxiso2(ib),nyiso2(ib),nx,ny,nx1,ny1,2)
        enddo
        boundary_filename = 
     .   '../'//boundary_filename(1:lnblnk(boundary_filename))//'.2'
        write(*,*) 'Writing ',trim(boundary_filename)
        open(35,file=boundary_filename)
        write(35,'(a)') '&BOUNDARY'
        write(35,'(1x,a,i2,a)') 'NBC=', nbc1,','
        line = 'BCCHAR= '
        len=8
        do ib = 1, nbc1
          line = line(1:len)//'  '''//bcchar(ib)//''','
          len=len+3+1+2
        end do
        write(35,'(1x,a)') line
        write(hlp_frm,'(a,i2,a)') '(1x,a,',nbc1,'(1x,i4,'',''))'
        write(35, hlp_frm) 'BCPOS=  ',(bcpos1(ib)  ,ib=1,nbc1)
        write(35, hlp_frm) 'BCSTART=',(bcstart1(ib),ib=1,nbc1)
        write(35, hlp_frm) 'BCEND=  ',(bcend1(ib)  ,ib=1,nbc1)
        write(35, hlp_frm) 'BCENE=  ',(bcene(ib)  ,ib=1,nbc1)
        write(35, hlp_frm) 'BCENI=  ',(bceni(ib)  ,ib=1,nbc1)
        write(35, hlp_frm) 'BCPOT=  ',(bcpot(ib)  ,ib=1,nbc1)
        write(hlp_frm,'(a,i2,a)') '(1x,a,1p,',nbc1,'(1x,g9.2,'',''))'
        write(35, hlp_frm) 'ENEPAR(1,1)=',(enepar(ib,1),ib=1,nbc1)
        if (maxval(abs(enepar(1:nbc1,2))).gt.0.0_R8) 
     >   write(35, hlp_frm) 'ENEPAR(1,2)=',(enepar(ib,2),ib=1,nbc1)
        write(35, hlp_frm) 'ENIPAR(1,1)=',(enipar(ib,1),ib=1,nbc1)
        if (maxval(abs(enipar(1:nbc1,2))).gt.0.0_R8) 
     >   write(35, hlp_frm) 'ENIPAR(1,2)=',(enipar(ib,2),ib=1,nbc1)
        write(35, hlp_frm) 'POTPAR(1,1)=',(potpar(ib,1),ib=1,nbc1)
        if (maxval(abs(potpar(1:nbc1,2))).gt.0.0_R8) 
     >   write(35, hlp_frm) 'POTPAR(1,2)=',(potpar(ib,2),ib=1,nbc1)
        write(hlp_frb,'(a,i2,a)') '(1x,a,i2,a,',ns1,'(1x,i3,'',''))'
        write(hlp_frp,'(a,i2,a)') 
     .   '(1x,a,i2,a,1p,',ns1,'(1x,g9.2,'',''))'
        do ib = 1, nbc1
          write(35, hlp_frb) 
     .     'BCCON(0,',ib,')=',(bccon(is,ib),is=0,ns1-1)
        end do
        do ib = 1, nbc1
          write(35, hlp_frp) 
     .     'CONPAR(0,',ib,',1)=',(conpar(is,ib,1),is=0,ns1-1)
          if (maxval(abs(conpar(0:ns1-1,ib,2))).gt.0.0_R8) 
     >     write(35, hlp_frp) 
     >    'CONPAR(0,',ib,',2)=',(conpar(is,ib,2),is=0,ns1-1)
          if (maxval(abs(conpar(0:ns1-1,ib,3))).gt.0.0_R8) 
     >     write(35, hlp_frp) 
     >    'CONPAR(0,',ib,',3)=',(conpar(is,ib,3),is=0,ns1-1)
        end do
        do ib = 1, nbc1
          write(35, hlp_frb) 
     .   'BCMOM(0,',ib,')=',(bcmom(is,ib),is=0,ns1-1)
        end do
        do ib = 1, nbc1
          write(35, hlp_frp) 
     .   'MOMPAR(0,',ib,',1)=',(mompar(is,ib,1),is=0,ns1-1)
          if (maxval(abs(mompar(0:ns1-1,ib,2))).gt.0.0_R8) 
     >     write(35, hlp_frp) 
     >    'MOMPAR(0,',ib,',2)=',(mompar(is,ib,2),is=0,ns1-1)
        end do
        write(35,'(1x,a,f10.7,a,f10.7,a)') 
     .   'GAMMAI=',gammai,', GAMMAE=',gammae,','
        write(35,'(1x,a,l1,a,l1,a)')
     .   'LBNDUSR=',lbndusr,', LFEEDBACK=',lfeedback,','
        write(35,'(1x,a,i2,a)') 'NNISO=',NNISO,','
        if (nniso.gt.0) then
          write(hlp_frm,'(a,i2,a)') '(1x,a,',nniso,'(1x,i2,'',''))'
          write(35, hlp_frm) 'NXISO1=',nxiso3(1:nniso)
          write(35, hlp_frm) 'NXISO2=',nxiso4(1:nniso)
          write(35, hlp_frm) 'NYISO1=',nyiso3(1:nniso)
          write(35, hlp_frm) 'NYISO2=',nyiso4(1:nniso)
          write(hlp_frm,'(a,i2,a)') '(1x,a,1p,',ns1,'(1x,g9.2,'',''))'
          write(35, hlp_frm) 'NIISO=',niiso(0:ns1-1)
          write(35, '(1x,a,1p,g9.2,a,g9.2,a,g9.2)')
     .     'TEISO=',teiso,', TIISO=',tiiso,', PHIISO=',phiiso
        endif
        if (boundary_time_mod.ne.0.0_R8) write(35,'(1x,a,f10.7,a)') 
     >   'BOUNDARY_TIME_MOD=', boundary_time_mod,','
        if (boundary_time_switch.ne.0.0_R8) write(35,'(1x,a,f10.7,a)') 
     >   'BOUNDARY_TIME_SWITCH=', boundary_time_switch,','
        if (maxval(lcbs).gt.0) then
          ncbs = 1
          do while (ncbs.lt.nbc1 .and. lcbs(ncbs+1).gt.0)
            ncbs = ncbs+1
          end do
          write(hlp_frm,'(a,i2,a)') '(1x,a,',ncbs,'(1x,i2,'',''))'
          write(35, hlp_frm) 'LCBS=',lcbs(1:ncbs)
        endif
        write(35,'(a)') '/'
        close(35)
      endif

* ..now read and convert the sputtering and wall data if they exist
      if(sputter_namelist_used.or.wall_namelist_used) then
*  ..find out all the data needed for initialisation
        call ipgetr ('b2news_BoRiS', BoRiS)
        call ipgeti ('b2stbr_sput_phys_col', sput_phys_col)
        call xertst (0.le.sput_phys_col.and.sput_phys_col.le.9,
     1   'faulty b2stbr_sput_phys_col')
        call ipgeti ('b2stbr_sput_phys_model', sput_phys_model)
        call ipgeti ('b2stbr_refl_model', refl_model)
        call ipgeti ('b2stbr_plate_model', plate_model)
        call ipgeti ('b2stbr_plate_option', plate_option)
        call xertst (1.le.plate_option .and.
     &               3.ge.plate_option, 'faulty b2stbr_plate_option')
        call ipgetr ('b2stbr_plate_temp', plate_temp)
        call xertst (0.0_R8.lt.plate_temp, 'faulty b2stbr_plate_temp')
        call ipgetr ('b2stbr_plate_thick', plate_thick)
        call xertst (0.0_R8.le.plate_thick, 'faulty b2stbr_plate_thick')
*  ..read in the namelist data
        if(sputter_namelist_used) then
          open(99,file='../b2.sputter_save.parameters')
          read(99,sputter_save)
          close(99)
*  ..convert the sputter namelist data
          sputter_yield_yt = 0.0_R8
          sputter_yield2_yt = 0.0_R8
          do ix1 = -1, nx1
            if (ix1.eq.-1) then
              ixa = -1
            else
              ixa = 2*nx-1-(4*nx1-1-2*ix1)*nx/(2*nx1)
            endif
            if (ix1.eq.nx1) then
              ixb = nx
            else
              ixb = (2*ix1+2*nx1+1)*nx/(2*nx1)-nx
            endif
            do iy1 = -1, ny1
              if (xymap1(ix1,iy1).gt.0) then
                if (iy1.eq.-1) then
                  iya = -1
                else
                  iya = 2*ny-1-(4*ny1-1-2*iy1)*ny/(2*ny1)
                endif
                if (iy1.eq.ny1) then
                  iyb = ny
                else
                  iyb = (2*iy1+2*ny1+1)*ny/(2*ny1)-ny
                endif
                do ix = ixa, ixb
                  do iy = iya, iyb
                    if (xymap(ix,iy).gt.0) then
                      sputter_tmp(0:ns-1)=sputter_yield(ix,iy,0:ns-1,1)
                      call sum_bundle
     .                 (ns1,ns,nsm,is_neutral1,bundle_frac,sputter_tmp)
                      sputter_yield_yt(ix1,iy1,0:ns1-1,1)=
     .                 sputter_yield_yt(ix1,iy1,0:ns1-1,1)+
     .                 sputter_tmp(0:ns1-1)
                      sputter_tmp(0:ns-1)=sputter_yield(ix,iy,0:ns-1,2)
                      call sum_bundle
     .                 (ns1,ns,nsm,is_neutral1,bundle_frac,sputter_tmp)
                      sputter_yield_yt(ix1,iy1,0:ns1-1,2)=
     .                 sputter_yield_yt(ix1,iy1,0:ns1-1,2)+
     .                 sputter_tmp(0:ns1-1)
                      sputter_tmp(0:ns-1)=sputter_yield2(ix,iy,0:ns-1,1)
                      call sum_bundle
     .                 (ns1,ns,nsm,is_neutral1,bundle_frac,sputter_tmp)
                      sputter_yield2_yt(ix1,iy1,0:ns1-1,1)=
     .                 sputter_yield2_yt(ix1,iy1,0:ns1-1,1)+
     .                 sputter_tmp(0:ns1-1)
                      sputter_tmp(0:ns-1)=sputter_yield2(ix,iy,0:ns-1,2)
                      call sum_bundle
     .                 (ns1,ns,nsm,is_neutral1,bundle_frac,sputter_tmp)
                      sputter_yield2_yt(ix1,iy1,0:ns1-1,2)=
     .                 sputter_yield2_yt(ix1,iy1,0:ns1-1,2)+
     .                 sputter_tmp(0:ns1-1)
                    end if
                  end do
                end do
              end if
            end do
          end do
*  ..write new sputtering namelist
          if(maxval(sputter_yield_yt).gt.0.0_R8.or.
     .       maxval(sputter_yield2_yt).gt.0.0_R8) then
            open(99,file='../b2.sputter_save.parameters.2')
            write(99,'(a13)') '&SPUTTER_SAVE'
            write(99,'(a8,i4,a1)') 'NX_SP = ', nx1, ','
            write(99,'(a8,i4,a1)') 'NY_SP = ', ny1, ','
            write(99,'(a8,i4,a1)') 'NS_SP = ', ns1, ','
            do ix = -1, nx1
              do iy = -1, ny1
                do is = 0, ns1-1
                  if ( sputter_yield_yt(ix,iy,is,1).ne.0.0_R8 .or.
     &                 sputter_yield_yt(ix,iy,is,2).ne.0.0_R8 .or.
     &                sputter_yield2_yt(ix,iy,is,1).ne.0.0_R8 .or.
     &                sputter_yield2_yt(ix,iy,is,2).ne.0.0_R8) then
                    do ib = 1, 2
                      write(99,'(a14,i4,a1,i4,a1,i3,a1,i1,a2,
     &                          1pe21.13e3,
     &                          a17,i4,a1,i4,a1,i3,a1,i1,a2,
     &                          1pe21.13e3,a1)')
     .                 'SPUTTER_YIELD(',ix,',',iy,',',is,',',ib,')=',
     .                  sputter_yield_yt(ix,iy,is,ib),
     .               ', SPUTTER_YIELD2(',ix,',',iy,',',is,',',ib,')=',
     .                  sputter_yield2_yt(ix,iy,is,ib),','
                    end do
                  end if
                end do
              end do
            end do
            write(99,'(a1)') '/'
            close(99)
          endif
        endif  ! sputter_namelist_used
*  ..convert some wall arrays
        if(wall_namelist_used) then
          call b2yt_wall(nx1,ny1,nx,ny,ns,ns1,is_neutral1,bundle_frac,
     .     plate_model,plate_temp,plate_thick)
*  ..write new wall namelist
          call b2yt_wall_final(nx1,ny1,ns1,
     .     plate_model,plate_temp,plate_thick)
        endif ! wall_namelist_used
      endif

*  ..write new b2.transport.parameters
      inquire(file='../b2.transport.parameters',exist=todo)
      if(.not.same_bundling.and.(transport_namelist.eq.1.or.todo)) then        
        call read_b2mod_transport_namelist
        write(*,*) 'Writing b2.transport.parameters.2'
        write(hlp_frm,'(a,i2,a)') 
     .   '(1x,a,i2,a,a,1p,',ns1,'(g9.2,'','',1x))'
        open(99,file='../b2.transport.parameters.2')
        write(99,'(a10)') '&TRANSPORT'
        if(flag_dna.ne.0) then
          cftmp(0:ns-1)=parm_dna(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          parm_dna(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 
     .     'FLAG_DNA=',flag_dna,',',' PARM_DNA=',parm_dna(0:ns1-1)
        endif
        if(flag_dpa.ne.0) then
          cftmp(0:ns-1)=parm_dpa(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          parm_dpa(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 
     .     'FLAG_DPA=',flag_dpa,',',' PARM_DPA=',parm_dpa(0:ns1-1)
        endif
        if(flag_vla.ne.0) then
          cftmp(0:ns-1)=parm_vla(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          parm_vla(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 
     .     'FLAG_VLA=',flag_vla,',',' PARM_VLA=',parm_vla(0:ns1-1)
        endif
        if(flag_vsa.ne.0) then
          cftmp(0:ns-1)=parm_vsa(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          parm_vsa(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 
     .     'FLAG_VSA=',flag_vsa,',',' PARM_VSA=',parm_vsa(0:ns1-1)
        endif
        if(flag_hci.ne.0) then
          cftmp(0:ns-1)=parm_hci(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          parm_hci(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 
     .     'FLAG_HCI=',flag_hci,',',' PARM_HCI=',parm_hci(0:ns1-1)
        endif
        if(flag_hce.ne.0) write(99,'(1x,a,i2,a,a,1p,1g9.2,a)') 
     .     'FLAG_HCE=',flag_hce,',',' PARM_HCE=',parm_hce,','
        if(flag_sig.ne.0) write(99,'(1x,a,i2,a,a,1p,1g9.2,a)') 
     .     'FLAG_SIG=',flag_sig,',',' PARM_SIG=',parm_sig,','
        if(flag_alf.ne.0) write(99,'(1x,a,i2,a,a,1p,1g9.2,a)') 
     .     'FLAG_ALF=',flag_alf,',',' PARM_ALF=',parm_alf,','
        if(transport_time_mod.gt.0.0_R8) then
          write(99,'(1x,a,a,a)')
     .     'TRANSPORT_FILENAME=''',trim(transport_filename),''','
          write(99,'(1x,a,1p,1g9.2,a)')
     .     'TRANSPORT_TIME_MOD=',transport_time_mod,','  
          write(99,'(1x,a,1p,1g9.2,a)')
     .     'TRANSPORT_TIME_SWITCH=',transport_time_switch,',' 
        endif 
        if(cflme.ge.0.0_R8) write(99,'(1x,a,1p,1g9.2,a)')
     .   'CFLME=',cflme,','
        if(cflmi.ge.0.0_R8) write(99,'(1x,a,1p,1g9.2,a)')
     .   'CFLMI=',cflmi,','
        if(cflmv.ge.0.0_R8) write(99,'(1x,a,1p,1g9.2,a)')
     .   'CFLMV=',cflmv,','
        write(99,'(a1)') '/'
        if(abs(maxval(vout_cnv(0:ns-1))).gt.0.0_R8) then
          write(hlp_frm,'(a,i2,a)') 
     .     '(1x,a,1p,',ns1,'(1x,g9.2,'',''))'
          cftmp(0:ns-1)=vout_cnv(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          vout_cnv(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 
     .     'VOUT_CNV=',vout_cnv(0:ns1-1)
          cftmp(0:ns-1)=pw0_cnv(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          pw0_cnv(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 
     .     'PW0_CNV=',pw0_cnv(0:ns1-1)
          cftmp(0:ns-1)=pw1_cnv(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          pw1_cnv(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 
     .     'PW1_CNV=',pw1_cnv(0:ns1-1)
          cftmp(0:ns-1)=pw2_cnv(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          pw2_cnv(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 
     .     'PW2_CNV=',pw2_cnv(0:ns1-1)
        endif
        close(99)
      endif
  
*  ..write new b2.transport.inputfile
      prdone = .false.
      allocate(px(1:ny1+2),pf(1:ny1+2))
      inquire(file='../b2.transport.inputfile',exist=todo)
      if(transport_inputfile.eq.1.or.todo) then
        ndata(:,:,:)=0
        tdata(:,:,:,:)=0.0_R8
        addspec(:,:,:)=-5
        region_flags=.true.
        no_pflux=.false.
        poloidal_scaling=.false.
        scaling_strength=0.0_R8
        scaling_width=0.0_R8
        scaling_ix_begin=-2
        scaling_ix_end=-2
        elm_time_begin=0.0_R8
        elm_time_end=0.0_R8
        elm_time_period=0.0_R8
        elm_ix_begin=-2
        elm_ix_end=-2
        transport_time_mod=0.0_R8
        transport_time_switch=0.0_R8
        allocate(n1data(1:nkind_data,1:nkind_coeff,0:nsm-1))
        allocate(t1data(3,1:DEF_NYD+2,1:nkind_coeff,0:nsm-1))
        allocate(ad1spec(ns1,1:nkind_coeff,0:ns1-1))
        allocate(pelm(1:ny1+2))
        n1data(:,:,:)=0
        t1data(:,:,:,:)=0.0_R8
        ad1spec(:,:,:)=-5
        open(99,file='../b2.transport.inputfile')
        read(99,transport)
        close(99)
        write(*,*) 'Writing b2.transport.inputfile.2'
        open(99,file='../b2.transport.inputfile.2')
        write(99,'(a10)') '&TRANSPORT'
        do j = 1, nkind_coeff
          do i = 1, nkind_data
cxpb  We do the mesh size transformation first
            do is = 0, ns-1
              ndat = ndata(i,j,is)
              if (ndat.gt.0 .and. is1(is).ge.0 .and. is1(is).lt.ns1) 
     .         then
                if (ndat.gt.ny1+2) then  ! need to coarsen the profile according to the new grid
                  if (.not.prdone) then
                    allocate(prdsa(1:ny1+2),pra(1:ny1+2))
                    allocate(r(1:ny+2),f(1:ny+2),felm(1:ny+2))
                    allocate(d(1:ny+2),delm(1:ny+2))
                    prdsa(1)=sqrt((cr1(jxa1,-1)-
     1               0.5_R8*(crx1(jxa1,-1,0)+crx1(jxa1,-1,1)))**2+
     2                (cz1(jxa1,-1)-
     3               0.5_R8*(cry1(jxa1,-1,0)+cry1(jxa1,-1,1)))**2)
                    do iy = 0, ny1
                      prdsa(iy+2)=prdsa(iy+1)+
     1                 sqrt((cr1(jxa1,iy)-cr1(jxa1,iy-1))**2+
     2                      (cz1(jxa1,iy)-cz1(jxa1,iy-1))**2)
                    enddo
                    ds_offset=(prdsa(jsep1+2)+prdsa(jsep1+3))/2.0_R8
                    do iy = 1, ny1+2
                      prdsa(iy) = prdsa(iy) - ds_offset
                    enddo 
                    prdone = .true.
                  endif
                  r(1:ndat) = tdata(1,1:ndat,j,is)
                  f(1:ndat) = tdata(2,1:ndat,j,is)
                  felm(1:ndat) = tdata(3,1:ndat,j,is)
                  call E01BEF(ndat,r,f,d,ifail)
                  call E01BEF(ndat,r,felm,delm,ifail)
                  nfitanf=1
                  do iy=1,ny1+2
                    pra(iy)=prdsa(iy)
                    if(pra(iy).lt.r(1)) nfitanf=nfitanf+1
                    nfitend=iy-1
                    if (pra(iy).GT.r(ndat)) exit
                    nfitend=iy
                  enddo
                  call xertst(nfitend.le.nfitanf, 
     .             'Could not complete transport profile!')
                  ndim=nfitend+1-nfitanf
                  do iy = 1, ndim
                    px(i) = pra(nfitanf-1+iy)
                  end do
                  call E01BFF(ndat,r,f,d,ndim,px,pf,ifail)
                  call E01BFF(ndat,r,felm,delm,ndim,px,pelm,ifail)
                else  ! simple case
                  ndim = ndat
                  px(1:ndim) = tdata(1,1:ndim,j,is)
                  pf(1:ndim) = tdata(2,1:ndim,j,is)
                  pelm(1:ndim) = tdata(3,1:ndim,j,is)
                endif
                n1data(i,j,is1(is))=ndim 
                t1data(1,1:ndim,j,is1(is))=px(1:ndim)   ! profile position
                t1data(2,1:ndim,j,is1(is))=pf(1:ndim)   ! profile value 
                t1data(3,1:ndim,j,is1(is))=pelm(1:ndim) ! ELM value
                write(99,'(1x,a,i1,a,i1,a,i2,a,i2,a)')
     .           'NDATA(',i,',',j,',',is1(is),')=',
     .            n1data(i,j,is1(is)),','
                iy = 0
cxpb Check if we need to include new additional species because of de-bundling
                if (bundle_frac(is1(is),is).lt.1.0_R8) then
                  do js = 0, ns1-1
                    if (js.ne.is1(is) .and. 
     &                  bundle_frac(js,is).ne.0.0_R8) then
                      iy = iy + 1
                      ad1spec(iy,j,is1(is)) = js
                      write(99,'(1x,a,i2,a,i1,a,i2,a,i2,a)') 
     .                 'ADDSPEC(',iy,',',j,',',is1(is),')=',
     .                  ad1spec(iy,j,is1(is)),','
                    endif
                  enddo
                endif
                do ix = 1, ns
cxpb Convert old additional species
                  if (addspec(ix,j,is).ge.0) then
                    todo = .true.
                    ib = 1
                    do while (ib.le.ns1 .and. 
     &                        ad1spec(ib,j,is1(is)).ne.-5)
                      todo = todo .and.
     &                 is1(addspec(ix,j,is)).ne.is1(is) .and.   ! avoid double-counting
     &                 is1(addspec(ix,j,is)).ne.ad1spec(ib,j,is1(is))
                      ib = ib + 1
                    enddo
                    if (todo) then
                      iy = iy + 1
                      ad1spec(iy,j,is1(is)) = is1(addspec(ix,j,is))
                      write(99,'(1x,a,i2,a,i1,a,i2,a,i2,a)') 
     .                 'ADDSPEC(',iy,',',j,',',is1(is),')=',
     .                  ad1spec(iy,j,is1(is)),','
                    endif
                  endif
                enddo
                do ix = 1, ndim
                  write(99,'(1x,a,i2,a,i1,a,i2,a,1p,3(g9.2,a))')
     .             'TDATA(1,',ix,',',j,',',is1(is),')=',
     .              t1data(1,ix,j,is1(is)),',',
     .              t1data(2,ix,j,is1(is)),',',
     .              t1data(3,ix,j,is1(is)),','
                enddo
              endif
            enddo
          enddo
          todo = .false.
          do ireg = 1, nnreg(0)
            todo = todo .or. .not.region_flags(ireg,j)
          enddo
          write(hlp_frm, '(a,i1,a)') '(1x,a,i1,a,',nnreg(0),'l1,a)'
          if (todo) write(99,hlp_frm)
     .     'REGION_FLAGS(1,',j,')=',
     .     (region_flags(ireg,j),ireg=1,nnreg(0)),','
        enddo
        write(99,'(1x,a,l1,a)') 'NO_PFLUX=',no_pflux,','
        if(transport_time_mod.gt.0.0_R8) then
          write(99,'(1x,a,a,a)')
     .     'TRANSPORT_FILENAME=''',trim(transport_filename),''','
          write(99,'(1x,a,1p,1g9.2,a)')
     .     'TRANSPORT_TIME_MOD=',transport_time_mod,','  
          write(99,'(1x,a,1p,1g9.2,a)')
     .     'TRANSPORT_TIME_SWITCH=',transport_time_switch,',' 
        endif 
        if (poloidal_scaling) write(99,'(1x,a,l1,a)') 
     .   'POLOIDAL_SCALING=',poloidal_scaling,','
        if (scaling_ix_begin.ne.-2) then
          scaling_ix_begin = conformal_convert(non_commensurate,1,
     .     scaling_ix_begin,jsep,nx,ny,nx1,ny1,1)
          write(99,'(1x,a,i4,a)')
     .     'SCALING_IX_BEGIN=',scaling_ix_begin,','
        endif
        if (scaling_ix_end.ne.-2) then
          scaling_ix_end = conformal_convert(non_commensurate,1,
     .     scaling_ix_end,jsep,nx,ny,nx1,ny1,3)
          write(99,'(1x,a,i4,a)')
     .     'SCALING_IX_END=',scaling_ix_end,','
        endif
        if (scaling_strength.ne.0.0_R8) write(99,'(1x,a,1p,g9.2,a)')
     .     'SCALING_STRENGTH=',scaling_strength,','
        if (scaling_width.ne.0.0_R8) write(99,'(1x,a,1p,g9.2,a)')
     .     'SCALING_WIDTH=',scaling_width,','
        if (elm_time_begin.ne.0.0_R8) write(99,'(1x,a,1p,g9.2,a)')
     .     'ELM_TIME_BEGIN=',elm_time_begin,','
        if (elm_time_end.ne.0.0_R8) write(99,'(1x,a,1p,g9.2,a)')
     .     'ELM_TIME_END=',elm_time_end,','
        if (elm_time_period.ne.0.0_R8) write(99,'(1x,a,1p,g9.2,a)')
     .     'ELM_TIME_PERIOD=',elm_time_period,','
        if (elm_ix_begin.ne.-2) then
          elm_ix_begin = conformal_convert(non_commensurate,1,
     .     elm_ix_begin,jsep,nx,ny,nx1,ny1,1)
          write(99,'(1x,a,i4,a)')
     .     'ELM_IX_BEGIN=',elm_ix_begin,','
        endif
        if (elm_ix_end.ne.-2) then
          elm_ix_end = conformal_convert(non_commensurate,1,
     .     elm_ix_end,jsep,nx,ny,nx1,ny1,3)
          write(99,'(1x,a,i4,a)')
     .     'ELM_IX_END=',elm_ix_end,','
        endif
        write(99,'(a1)') '/'
        close(99)
        deallocate(t1data)
        deallocate(n1data) 
        deallocate(ad1spec)       
        if (prdone) deallocate(felm,delm)
        deallocate(pelm)
      endif

*  ..write new b2.transport_models_save.parameters
      inquire(file='../b2.transport_models_save.parameters',exist=todo)
      if (todo .and. (ny.ne.ny1 .or. non_commensurate.eq.1)) then
        call alloc_b2mod_transport_models(nx,ny,nsm)
        call read_b2mod_transport_models
        allocate(dsa(-1:ny), deta(-1:ny), dsa1(-1:ny1))
        allocate(eta_hce_mult1(-1:ny1))
        eta_hce_mult1 = 1.0_R8
        ifail = 0
        dsa(-1)=sqrt((cr(jxa,-1)-
     1                0.5_R8*(crx(jxa,-1,0)+crx(jxa,-1,1)))**2+
     2               (cz(jxa,-1)-
     3                0.5_R8*(cry(jxa,-1,0)+cry(jxa,-1,1)))**2)
        do iy = 0, ny
          dsa(iy)=dsa(iy-1)+
     1            sqrt((cr(jxa,iy)-cr(jxa,iy-1))**2+
     2                 (cz(jxa,iy)-cz(jxa,iy-1))**2)
        enddo
        ds_offset=(dsa(jsep)+dsa(jsep+1))/2.0_R8
        do iy = -1, ny
          dsa(iy) = dsa(iy) - ds_offset
        enddo 
        dsa1(-1)=sqrt((cr1(jxa1,-1)-
     1                 0.5_R8*(crx1(jxa1,-1,0)+crx1(jxa1,-1,1)))**2+
     2                (cz1(jxa1,-1)-
     3                 0.5_R8*(cry1(jxa1,-1,0)+cry1(jxa1,-1,1)))**2)
        do iy = 0, ny1
          dsa1(iy)=dsa1(iy-1)+
     1             sqrt((cr1(jxa1,iy)-cr1(jxa1,iy-1))**2+
     2                  (cz1(jxa1,iy)-cz1(jxa1,iy-1))**2)
        enddo
        ds_offset=(dsa1(jsep1)+dsa1(jsep1+1))/2.0_R8
        do iy = -1, ny1
          dsa1(iy) = dsa1(iy) - ds_offset
        enddo 
        call E01BEF(ny+2,dsa,eta_hce_mult,deta,ifail)
        call E01BFF(ny+2,dsa,eta_hce_mult,deta,
     .              ny1+2,dsa1,eta_hce_mult1,ifail)
        if (maxval(eta_hce_mult1).ne.1.0_R8 .or. 
     .      minval(eta_hce_mult1).ne.1.0_R8) then
          write(*,*) 'Writing b2.transport_models_save.parameters.2'
          open(99,file='../b2.transport_models_save.parameters.2')
          write(99,'(a9)') '&TRANSPORT_MODELS'
          write(hlp_frm,'(a,i4,a)') '(1x,a,1p,',ny1+2,'(1x,g9.2,'',''))'
          write(99,hlp_frm) 
     .     'ETA_HCE_MULT=',eta_hce_mult1(-1:ny1)
          write(99,'(a1)') '/'
          close(99)
        endif
        deallocate(dsa,dsa1,deta)
        deallocate(eta_hce_mult1)
      endif

*  ..write new b2.sources.profile
      inquire(file='../b2.sources.profile',exist=todo)
      if(sources_inputfile.eq.1.or.todo) then
        nsdata(:,:,:)=0
        sdata(:,:,:,:)=0.0_R8
        divheat=0.0_R8
        allocate(ns1data(1:nkind_data,1:nkind_source,0:nsm-1))
        allocate(s1data(2,1:DEF_NYD+2,1:nkind_source,0:nsm-1))
        ns1data(:,:,:)=0
        s1data(:,:,:,:)=0.0_R8
        open (92,file='../b2.sources.profile')
        read(92,profile)
        close(92)
        write(*,*) 'Writing b2.sources.profile.2'
        open(99,file='../b2.sources.profile.2')
        write(99,'(a9)') '&PROFILE'
        do i = 1, nkind_data 
          do j = 1, nkind_source
            do is = 0, ns-1
              ndat=nsdata(i,j,is)
              if (ndat.gt.0.and.is1(is).ge.0.and.is1(is).lt.ns1) then
                clash(1:ns-1) = .false.
                shared(0:ns-1)= .false.
                w1d(0)=bundle_frac(is1(is),is)
                ks = 0
                do js = 0, ns-1
                  if (js.ne.is .and. nsdata(i,j,js).gt.0 .and.
     &                bundle_frac(is1(is),js).ne.0.0_R8) then  ! must check the profiles do not clash
                    shared(js)=.true.
                    ks = ks + 1
                    w1d(ks)=bundle_frac(is1(is),js)
cxpb Do we have the same number of points?
                    nds(ks) = nsdata(i,j,js)
                    clash(ks) = clash(ks) .or. ndat.ne.nds(ks)
                    if (.not.clash(ks)) then
cxpb Are the points in the same place?
                      do ix = 1, ndat
                        clash(ks) = clash(ks) .or.
     &                   abs(sdata(1,ix,j,is)-sdata(1,ix,j,js)).le.
     .                   1.0e-6_R8
                      enddo
                    endif
                  endif
                enddo
                ik = 0
                inc = 0
                s1data(1,1:ndat,j,is)=sdata(1,1:ndat,j,is)
                s1data(2,1:ndat,j,is)=sdata(2,1:ndat,j,is)*w1d(0)
                do js = 0, ns-1
                  if (shared(js)) then
                    ik = ik+1
                    if (.not.clash(ik)) then
cxpb We can sum up all the profiles that do not clash
                      inc = inc+1
                      s1data(2,1:ndat,j,is)=sdata(2,1:ndat,j,is)+
     +                                      sdata(2,1:ndat,j,js)*w1d(ik)
                    endif
                  endif
                enddo
cxpb Some profiles clash: we must compute a new combined profile 
                if (inc.lt.ks) then
                  nnc = ks-inc
                  ik = 0
                  allocate(sr1(1:(ny+2)),sf1(1:(ny+2)))
                  allocate(sd1(1:(ny+2)))
                  sr1(1:ndat) = sdata(1,1:ndat,j,is)
                  sf1(1:ndat) = sdata(2,1:ndat,j,is)
                  call E01BEF(ndat,sr1,sf1,sd1,ifail)
cxpb Calculate radial locations of combined profile spr(1:ndim)
                  allocate(spr(1:nnc*(ny+2)))
                  ndim = ndat
                  spr(1:ndat) = sr1(1:ndat)
                  do js = 0, ns-1
                    if (shared(js)) then 
                      ik = ik+1
                      if (clash(ik)) then
                        i1 = 1
                        i2 = 1
                        do while(i2.le.nds(ik)) 
                          if (abs(spr(i1)-sr2(i2)).lt.1.0e-6_R8) then  ! shared point
                            i1 = i1+1
                            i2 = i2+1
                          else if (spr(i1).lt.sr2(i2)) then  ! catch up with 2nd profile
                            i1 = i1+1
                          else ! add a new point to the profile
                            do iy=ndim,i1,-1
                              spr(iy+1)=spr(iy)
                            enddo
                            spr(i1)=sr2(i2)
                            ndim=ndim+1
                            i2 = i2+1  
                          endif
                        enddo
                      endif
                    endif
                  enddo
                  allocate(spf1(1:ndim))
                  nfitanf=1
                  do iy = 1, ndim
                    if(spr(iy).lt.sr1(1)) nfitanf=nfitanf+1
                    nfitend=iy-1
                    if(spr(iy).gt.sr1(ndat)) exit
                    nfitend=iy
                  enddo
                  call xertst(nfitend.le.nfitanf, 
     .             'Could not complete initial source profile!')
                  ndm2=nfitend+1-nfitanf
                  call E01BFF(ndat,sr1,sf1,sd1,ndm2,
     .             spr(nfitanf:nfitend),spf1(nfitanf:nfitend),ifail)
                  spf1(1:nfitanf-1)=0.0_R8
                  spf1(nfitend+1:ndim)=0.0_R8
                  allocate(s1d(2,1:ndim,1:nkind_source,0:nsm-1))
                  s1d(1,1:ndim,j,is)=spr(1:ndim)
                  s1d(2,1:ndim,j,is)=spf1(1:ndim)*w1d(0)
                  deallocate(sr1,sf1,sd1,spf1)
                  nsdata(i,j,is)=ndim
                  ndat=ndim
                  ik = 0
                  do js = 0, ns-1  
                    if (shared(js)) then 
                      ik = ik+1
                      if (clash(ik)) then
                        allocate(sr2(1:nds(ik)),sf2(1:nds(ik)))
                        allocate(sd2(1:nds(ik)))
                        allocate(spf2(1:ndim))
                        sr2(1:nds(ik)) = s1data(1,1:nds(ik),j,js)
                        sf2(1:nds(ik)) = s1data(2,1:nds(ik),j,js)
                        nfitanf=1
                        do iy = 1, ndim
                          if(spr(iy).lt.sr2(1)) nfitanf=nfitanf+1
                          nfitend=iy-1
                          if(spr(iy).gt.sr2(nds(ik))) exit
                          nfitend=iy
                        enddo
                        call xertst(nfitend.le.nfitanf, 
     .                'Could not complete additional source profile!')
                        ndm2=nfitend+1-nfitanf
                        call E01BEF(nds(ik),sr2,sf2,sd2,ifail)
                        call E01BFF(nds(ik),sr2,sf2,sd2,ndm2,
     .                   spr(nfitanf:nfitend),spf2(nfitanf:nfitend),
     .                   ifail)
                        spf2(1:nfitanf-1)=0.0_R8
                        spf2(nfitend+1:ndim)=0.0_R8
                        s1d(2,1:ndim,j,is)=
     .                   s1d(2,1:ndim,j,is)+spf2(1:ndim)*w1d(ik)
                        deallocate(sr2,sf2,sd2,spf2)
                      endif
                    endif
                  enddo
                  deallocate(spr)
                else  ! no clash : copy over data
                  allocate(s1d(2,1:ndat,1:nkind_source,0:nsm-1))
                  s1d(1,1:ndat,j,is)=s1data(1,1:ndat,j,is)
                  s1d(2,1:ndat,j,is)=s1data(2,1:ndat,j,is)
                endif
                if (bundle_frac(is1(is),js).eq.1.0_R8) 
     >            nsdata(i,j,js)=0   ! This profile has been absorbed in the current one
                if (ndat.gt.ny1+2) then  ! need to coarsen the profile according to the new grid
                  if (.not.prdone) then
                    allocate(prdsa(1:ny1+2),pra(1:ny1+2))
                    allocate(r(1:ny+2),f(1:ny+2))
                    allocate(d(1:ny+2))
                    allocate(px(1:ny1+2),pf(1:ny1+2))
                    prdsa(1)=sqrt((cr1(jxa1,-1)-
     1               0.5_R8*(crx1(jxa1,-1,0)+crx1(jxa1,-1,1)))**2+
     2                (cz1(jxa1,-1)-
     3               0.5_R8*(cry1(jxa1,-1,0)+cry1(jxa1,-1,1)))**2)
                    do iy = 0, ny1
                      prdsa(iy+2)=prdsa(iy+1)+
     1                 sqrt((cr1(jxa1,iy)-cr1(jxa1,iy-1))**2+
     2                      (cz1(jxa1,iy)-cz1(jxa1,iy-1))**2)
                    enddo
                    ds_offset=(prdsa(jsep1+2)+prdsa(jsep1+3))/2.0_R8
                    do iy = 1, ny1+2
                      prdsa(iy) = prdsa(iy) - ds_offset
                    enddo 
                    prdone = .true.
                  endif
                  r(1:ndat) = s1d(1,1:ndat,j,is)
                  f(1:ndat) = s1d(2,1:ndat,j,is)
                  call E01BEF(ndat,r,f,d,ifail)
                  nfitanf=1
                  do iy=1,ny1+2
                    pra(iy)=prdsa(iy)
                    if(pra(iy).lt.r(1)) nfitanf=nfitanf+1
                    nfitend=iy-1
                    if (pra(iy).GT.r(ndat)) exit
                    nfitend=iy
                  enddo
                  call xertst(nfitend.le.nfitanf, 
     .             'Could not complete source profile!')
                  ndim=nfitend+1-nfitanf
                  do iy = 1, ndim
                    px(iy) = pra(nfitanf-1+iy)
                  end do
                  call E01BFF(ndat,r,f,d,ndim,px,pf,ifail)
                else  ! simple case
                  ndim = ndat
                  px(1:ndim) = s1d(1,1:ndim,j,is)
                  pf(1:ndim) = s1d(2,1:ndim,j,is)
                endif
                ns1data(i,j,is1(is))=ndim
                s1data(1,1:ndim,j,is1(is))=px(1:ndim)  ! profile position
                s1data(2,1:ndim,j,is1(is))=pf(1:ndim)  ! profile value
                write(99,'(1x,a,i1,a,i1,a,i2,a,i2,a)')
     .           'NSDATA(',i,',',j,',',is1(is),')=',
     .            ns1data(i,j,is1(is)),','
                do ix = 1, ndim
                  write(99,'(1x,a,i2,a,i1,a,i2,a,1p,2(g9.2,a))')
     .             'SDATA(1,',ix,',',j,',',is1(is),')=',
     .             s1data(1,ix,j,is1(is)),',',s1data(2,ix,j,is1(is)),','
                enddo
                write(99,1028)
              endif
            enddo
          enddo
        enddo
        if (divheat.ne.0.0_R8) write(99,'(1x,a,1p,1g9.2,a)') 
     .   'DIVHEAT=',divheat,','
        write(99,'(a1)') '/'
        close(99)
        deallocate(s1data,s1d)
        deallocate(ns1data)
      endif
      if (prdone) deallocate(prdsa,pra,r,f,d)
      deallocate(px,pf)
     
*  ..write new b2.numerics.parameters
      inquire(file='../b2.numerics.parameters',exist=todo)
      if(.not.same_bundling.and.(numerics_namelist.ge.1.or.todo)) then
        call read_b2mod_numerics_namelist
        write(*,*) 'Writing b2.numerics.parameters.2'
        open(99,file='../b2.numerics.parameters.2')
        write(99,'(a9)') '&NUMERICS'
        write(hlp_frm,'(a,i2,a)') 
     .   '(1x,a,i2,a,1p,',ns1,'(1x,g9.2,'',''))'
        if(minval(dtco(0:ns-1,0:nnreg(0))).ne.1.0_R8 .or.
     &     maxval(dtco(0:ns-1,0:nnreg(0))).ne.1.0_R8) then
          do ireg = 0, nnreg(0)
            cftmp(0:ns-1)=dtco(0:ns-1,ireg)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            dtco(0:ns1-1,ireg)=cftmp(0:ns1-1)
            write(99,hlp_frm)
     .       'DTCO(0,',ireg,')=',(dtco(is,ireg),is=0,ns1-1)
          enddo
        endif
        if(minval(dtmo(0:ns-1,0:nnreg(0))).ne.1.0_R8 .or.
     &     maxval(dtmo(0:ns-1,0:nnreg(0))).ne.1.0_R8) then
          do ireg = 0, nnreg(0)
            cftmp(0:ns-1)=dtmo(0:ns-1,ireg)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            dtmo(0:ns1-1,ireg)=cftmp(0:ns1-1)
            write(99,hlp_frm)
     .       'DTMO(0,',ireg,')=',(dtmo(is,ireg),is=0,ns1-1)
          enddo
        endif
        write(hlp_frm,'(a,i2,a)') 
     .   '(1x,a,1p,',nnreg(0)+1,'(1x,g9.2,'',''))'
        if(minval(dtee(0:nnreg(0))).ne.1.0_R8 .or.
     &     maxval(dtee(0:nnreg(0))).ne.1.0_R8) write(99,hlp_frm)
     .   'DTEE=',(dtee(ireg),ireg=0,nnreg(0))
        if(minval(dtei(0:nnreg(0))).ne.1.0_R8 .or.
     &     maxval(dtei(0:nnreg(0))).ne.1.0_R8) write(99,hlp_frm)
     .   'DTEI=',(dtei(ireg),ireg=0,nnreg(0))
        write(hlp_frm,'(a,i2,a)') '(1x,a,i2,a,',ns1,'l1,a)'
        do ireg = 0, nnreg(0)
          todo = .false.
          do is = 0, ns-1
            todo = todo .or. .not.solveco(is,ireg)
          enddo
          if (todo) then
            call convert_bundle_l(ns1,ns,nsm,is_neutral1,bundle_frac,
     .       solveco(0,ireg))
            write(99,hlp_frm)
     .       'SOLVECO(0,',ireg,')=',(solveco(is,ireg),is=0,ns1-1)
          endif
        enddo
        do ireg = 0, nnreg(0)
          todo = .false.
          do is = 0, ns-1
            todo = todo .or. .not.solvemo(is,ireg)
          enddo
          if (todo) then
            call convert_bundle_l(ns1,ns,nsm,is_neutral1,bundle_frac,
     .       solvemo(0,ireg))
            write(99,hlp_frm)
     .       'SOLVEMO(0,',ireg,')=',(solvemo(is,ireg),is=0,ns1-1)
          endif
        enddo
        write(hlp_frm,'(a,i2,a)') '(1x,a,',nnreg(0)+1,'l1,a)'
        todo = .false.
        do ireg = 0, nnreg(0)
          todo = todo .or. .not.solvepo(ireg)
        enddo
        if (todo) write(99,hlp_frm)
     .    'SOLVEPO=',(solvepo(ireg),is=0,nnreg(0))
        todo = .false.
        do ireg = 0, nnreg(0)
          todo = todo .or. .not.solveee(ireg)
        enddo
        if (todo) write(99,hlp_frm)
     .    'SOLVEEE=',(solveee(ireg),is=0,nnreg(0))
        todo = .false.
        do ireg = 0, nnreg(0)
          todo = todo .or. .not.solveei(ireg)
        enddo
        if (todo) write(99,hlp_frm)
     .    'SOLVEEI=',(solveei(ireg),is=0,nnreg(0))
        if(numerics_time_mod.gt.0.0_R8) then
          write(99,'(1x,a,a,a)')
     .     'NUMERICS_FILENAME=''',trim(numerics_filename),''','
          write(99,'(1x,a,1p,1g9.2,a)')
     .     'NUMERICS_TIME_MOD=',numerics_time_mod,','  
          write(99,'(1x,a,1p,1g9.2,a)')
     .     'NUMERICS_TIME_SWITCH=',numerics_time_switch,',' 
        endif 
        if(core_dt_factor.ne.1.0_R8) write(99,'(1x,a,1p,1g9.2,a)')
     .   'CORE_DT_FACTOR=',core_dt_factor,','
        if(core_dt_suppression.ne.1.0_R8) write(99,'(1x,a,1p,1g9.2,a)')
     .   'CORE_DT_SUPPRESSION=',core_dt_suppression,','
        if(time_factor_required.ne.0.1_R8) write(99,'(1x,a,1p,1g9.2,a)')
     .   'TIME_FACTOR_REQUIRED=',time_factor_required,','
        write(99,'(a1)') '/'
        close(99)
      endif

*  ..write new b2.atomic_physics_rescale.parameters
      inquire(file='../b2.atomic_physics_rescale.parameters',exist=todo)
      if(atomic_physics_rescale_flag.ne.0.or.todo) then
        rescale_sa=1.0_R8
        rescale_ra=1.0_R8
        rescale_qa=1.0_R8
        rescale_cx=1.0_R8
        rescale_rd=1.0_R8
        rescale_br=1.0_R8
        open(99,file='../b2.atomic_physics_rescale.parameters')
        write(*,*) 'Reading ../b2.atomic_physics_rescale.parameters'
        read(99,atomic_physics_rescale)
        write(*,atomic_physics_rescale)
        close(99)
        write(*,*) 'Read ../b2.atomic_physics_rescale.parameters'
        call xertst (min(
     1   minval(rescale_sa(0:ns-1)),minval(rescale_ra(0:ns-1)),
     2   minval(rescale_qa(0:ns-1)),minval(rescale_cx(0:ns-1)),
     3   minval(rescale_rd(0:ns-1)),minval(rescale_br(0:ns-1)))
     4   .gt.0.0_R8,
     4   'faulty value set in b2.atomic_physics_rescale.parameters')
        write(*,*) 'Writing b2.atomic_physics_rescale.parameters.2'
        write(hlp_frm,'(a,i2,a)') '(1x,a,1p,',ns1,'(1x,g9.2,'',''))'
        open(99,file='../b2.atomic_physics_rescale.parameters.2')
        write(99,'(a10)') '&ATOMIC_PHYSICS_RESCALE'
        if (maxval(rescale_sa(0:ns-1)).ne.1.0_R8 .or.
     &      minval(rescale_sa(0:ns-1)).ne.1.0_R8) then
          cftmp(0:ns-1)=rescale_sa(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          rescale_sa(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 'RESCALE_SA=',rescale_sa(0:ns1-1)
        endif        
        if (maxval(rescale_ra(0:ns-1)).ne.1.0_R8 .or.
     &      minval(rescale_ra(0:ns-1)).ne.1.0_R8) then
          cftmp(0:ns-1)=rescale_ra(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          rescale_ra(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 'RESCALE_RA=',rescale_ra(0:ns1-1)
        endif        
        if (maxval(rescale_qa(0:ns-1)).ne.1.0_R8 .or.
     &      minval(rescale_qa(0:ns-1)).ne.1.0_R8) then
          cftmp(0:ns-1)=rescale_qa(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          rescale_qa(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 'RESCALE_QA=',rescale_qa(0:ns1-1)
        endif        
        if (maxval(rescale_cx(0:ns-1)).ne.1.0_R8 .or.
     &      minval(rescale_cx(0:ns-1)).ne.1.0_R8) then
          cftmp(0:ns-1)=rescale_cx(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          rescale_cx(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 'RESCALE_CX=',rescale_cx(0:ns1-1)
        endif        
        if (maxval(rescale_rd(0:ns-1)).ne.1.0_R8 .or.
     &      minval(rescale_rd(0:ns-1)).ne.1.0_R8) then
          cftmp(0:ns-1)=rescale_rd(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          rescale_rd(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 'RESCALE_RD=',rescale_rd(0:ns1-1)
        endif        
        if (maxval(rescale_br(0:ns-1)).ne.1.0_R8 .or.
     &      minval(rescale_br(0:ns-1)).ne.1.0_R8) then
          cftmp(0:ns-1)=rescale_br(0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          rescale_br(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm) 'RESCALE_BR=',rescale_br(0:ns1-1)
        endif        
        close(99)
      endif

*  ..write new input.dat file if present
      IUNIN = 86
      IUNOT = 87
      if (open_file(IUNIN,'input.dat','')) then
        read(IUNIN,'(A72)',end=8687) ZEILE   ! TXTRUN 
        open(unit=IUNOT,file='input.dat.2',status='unknown')
        write(*,*) 'Writing input.dat.2'
        write(IUNOT,'(A72)') ZEILE
*** 1. DATA FOR OPERATING MODE
  109   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1).EQ.'*') THEN
          WRITE(IUNOT,'(A72)') ZEILE
          GOTO 109
        ELSE
#ifdef NEW_SIZE
          READ (ZEILE,6666) idum(0:7)  ! NMACH,NMODE,NTCPU,NFILE,NITER0,NITER,
                                       ! NTIME0,NTIME
          NMODE = idum(1)
#else
          READ (ZEILE,6666) idum(0:5)  ! NMACH,NMODE,NTIME,NFILE,NITER,lfrm15
          NMODE = idum(1)
#endif
          WRITE(IUNOT,'(A72)') ZEILE
        ENDIF
        READ (IUNIN,6665) ldum(0:4)  ! NLSCL,NLTEST,NLANA,NLDRFT,NLCRR
        WRITE(IUNOT,6665) ldum(0:4)
*** 2. DATA FOR STANDARD MESH
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
  201   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 201
        READ (ZEILE,6666) (INDGRD(J),J=1,3)
C  RADIAL MESH
        READ (IUNIN,6665) ldum(0)    ! nlrad
        WRITE(IUNOT,6665) ldum(0)
        IF (ldum(0)) THEN
          READ (IUNIN,6665) ldum(0:6)  ! NLSLB,NLCRC,NLELL,NLTRI,NLPLG,NLFEM,NLGEN
          WRITE(IUNOT,6665) ldum(0:6)
#ifdef NEW_SIZE
          READ (IUNIN,6666) idum(0:4)  ! NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT
          NR1ST = idum(0)
          NPPLG = idum(3)
          if (ldum(0).or.ldum(1).or.ldum(2).or.ldum(3).or.ldum(4)) then
            idum(0) = ny1+1
            IF (INDGRD(1).LE.4) idum(1) = jsep1+1
          endif
          idum(5:6) = 0
#else
          READ (IUNIN,6666) idum(0:6)  ! NR1ST,NRSEP,NRPLG,NPPLG,NRFOU,NRKNOT,NRTRI
          NR1ST = idum(0)
          NPPLG = idum(3)
          if (ldum(0).or.ldum(1).or.ldum(2).or.ldum(3).or.ldum(4)) then
            idum(0) = ny1+1
            IF (INDGRD(1).LE.4) idum(1) = jsep1+1
          endif
#endif
          WRITE(IUNOT,6666) idum(0:6)
          IF (INDGRD(1).LE.5) THEN
            IF (ldum(0).OR.ldum(1).OR.ldum(2).OR.ldum(3)) THEN
              READ (IUNIN,'(A72)') ZEILE ! RIA,RGA,RAA,RRA
              WRITE(IUNOT,'(A72)') ZEILE 
              IF (ldum(2).OR.ldum(3)) THEN
                READ (IUNIN,'(A72)') ZEILE  ! EP1IN,EP1OT,EP1CH,EXEP1
                WRITE(IUNOT,'(A72)') ZEILE
                READ (IUNIN,'(A72)') ZEILE  ! ELLIN,ELLOT,ELLCH,EXELL
                WRITE(IUNOT,'(A72)') ZEILE
                IF (ldum(3)) THEN
                  READ (IUNIN,'(A72)') ZEILE  ! TRIIN,TRIOT,TRICH,EXTRI
                  WRITE(IUNOT,'(A72)') ZEILE
                ENDIF
              ENDIF
            ENDIF
            IF (ldum(4)) THEN
              READ (IUNIN,'(A72)') ZEILE  ! XPCOR,YPCOR,ZPCOR,PLREFL
              WRITE(IUNOT,'(A72)') ZEILE
              allocate(NPOINT(2,idum(3)))
              allocate(XPOL(NR1ST,idum(2)),YPOL(NR1ST,idum(2)))
#ifdef NEW_SIZE
              READ (IUNIN,6666) (NPOINT(1,K),NPOINT(2,K),K=1,idum(3))
              WRITE(IUNOT,6666) (NPOINT(1,K),NPOINT(2,K),K=1,idum(3))
              DO 212 I=1,NR1ST
                READ (IUNIN,6664) (XPOL(I,J),YPOL(I,J),J=1,idum(2))
                WRITE(IUNOT,6664) (XPOL(I,J),YPOL(I,J),J=1,idum(2))
212           CONTINUE
#else
              DO 211 I=1,NR1ST
                READ (IUNIN,6666) (NPOINT(1,K),NPOINT(2,K),K=1,idum(3))
                WRITE(IUNOT,6666) (NPOINT(1,K),NPOINT(2,K),K=1,idum(3))
                READ (IUNIN,6664) (XPOL(I,J),YPOL(I,J),    J=1,idum(2))
                WRITE(IUNOT,6664) (XPOL(I,J),YPOL(I,J),    J=1,idum(2))
211           CONTINUE
#endif
              deallocate(NPOINT)
              deallocate(XPOL,YPOL)
            ENDIF
#ifdef NEW_SIZE
            IF (ldum(5)) THEN  ! NLFEM
              READ (IUNIN,6666) idum(0)  ! NRKNOT
              WRITE(IUNOT,6666) idum(0)
              allocate(XTRIAN(idum(0)),YTRIAN(idum(0)))
              READ (IUNIN,6664) (XTRIAN(I),I=1,idum(0))
              WRITE(IUNOT,6664) (XTRIAN(I),I=1,idum(0))
              READ (IUNIN,6664) (YTRIAN(I),I=1,idum(0))
              WRITE(IUNOT,6664) (YTRIAN(I),I=1,idum(0))
              deallocate(XTRIAN,YTRIAN)
              DO 213 I=1,NR1ST
                READ (IUNIN,6666) idum(0:3) ! J,NECKE (1,J),NECKE (2,J),NECKE(3,J)
                WRITE(IUNOT,6666) idum(0:3) ! J,NECKE (1,J),NECKE (2,J),NECKE(3,J)
                READ (IUNIN,6666) idum(1:3) !   NCHBAR(1,J),NSEITE(1,J),INMTI(1,J)
                WRITE(IUNOT,6666) idum(1:3) !   NCHBAR(1,J),NSEITE(1,J),INMTI(1,J)
                READ (IUNIN,6666) idum(1:3) !   NCHBAR(2,J),NSEITE(2,J),INMTI(2,J)
                WRITE(IUNOT,6666) idum(1:3) !   NCHBAR(2,J),NSEITE(2,J),INMTI(2,J)
                READ (IUNIN,6666) idum(1:3) !   NCHBAR(3,J),NSEITE(3,J),INMTI(3,J)
                WRITE(IUNOT,6666) idum(1:3) !   NCHBAR(3,J),NSEITE(3,J),INMTI(3,J)
213           CONTINUE
            ENDIF
#endif
          ELSEIF (INDGRD(1).EQ.6) THEN
            READ (IUNIN,'(A72)') ZEILE
            WRITE(IUNOT,'(A72)') ZEILE
          ENDIF
        ENDIF
C  POLOIDAL MESH
  220   READ (IUNIN,'(A72)') ZEILE  ! NLPOL
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 220
#ifndef NEW_SIZE
        READ (IUNIN,6665) ldum(0:3) ! NLPLY,NLPLA,NLPLP,NLPLF
        WRITE(IUNOT,6665) ldum(0:3)
        READ (IUNIN,6666) idum(0:1) ! NP2ND,NPSEP
        NP2ND = idum(0)
        idum(0) = nx1 + max(1,nncut*3)
        WRITE(IUNOT,6666) idum(0:1)
        IF (INDGRD(2).LE.5) THEN
          READ (IUNIN,'(A72)') ZEILE ! YIA,YAGRID,YAA
          WRITE(IUNOT,'(A72)') ZEILE
        ELSEIF (INDGRD(2).EQ.6) THEN
        ENDIF
#else
        READ (IUNIN,6665) ldum(0:2) ! NLPLY,NLPLA,NLPLP
        WRITE(IUNOT,6665) ldum(0:2), .false.
        READ (IUNIN,6666) idum(0:3) ! NP2ND,NPSEP,NPPLA,NPPER
        NP2ND = idum(0)
        idum(0) = nx1 + max(1,nncut*3)
        WRITE(IUNOT,6666) idum(0:3)
        IF (INDGRD(2).LE.5) THEN
          READ (IUNIN,'(A72)') ZEILE ! YIA,YGA,YAA,YYA
          WRITE(IUNOT,'(A72)') ZEILE
        ELSEIF (INDGRD(2).EQ.6) THEN
        ENDIF
#endif
C  TOROIDAL MESH
  230   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 230
#ifndef NEW_SIZE
        READ (IUNIN,6665) ldum(0:3) ! NLTRZ,NLTRA,NLTRT,NLTRF
        WRITE(IUNOT,6665) ldum(0:3)
#else
        READ (IUNIN,6665) ldum(0:2) ! NLTRZ,NLTRA,NLTRT
        WRITE(IUNOT,6665) ldum(0:2), .false.
#endif        
        READ (IUNIN,6666) idum(0:3) ! NT3RD,NTSEP,NTTRA,NTPER
        WRITE(IUNOT,6666) idum(0:3)
        IF (INDGRD(3).LE.5) THEN
          READ (IUNIN,'(A72)') ZEILE ! ZIA,ZGA,ZAA,ZZA,ROA
          WRITE(IUNOT,'(A72)') ZEILE
        ELSEIF (INDGRD(3).EQ.6) THEN
        ENDIF
C  MESH MULTIPLICATION
  240   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 240
        READ (ZEILE,6665) ldum(0)
        IREAD=0
#ifndef NEW_SIZE
        READ (IUNIN,6666) idum(0) ! NBMLT
        WRITE(IUNOT,6666) idum(0)
        allocate(VOLCOR(idum(0)))
        IF (idum(0).GT.1) THEN
          READ (IUNIN,6664) (VOLCOR(I),I=1,idum(0))
          WRITE(IUNOT,6664) (VOLCOR(I),I=1,idum(0))
        ENDIF
        deallocate(VOLCOR)
#else
        IF (ldum(0)) THEN
          READ (IUNIN,6666) idum(0)
          WRITE(IUNOT,6666) idum(0)
          allocate(VOLCOR(idum(0)))
          READ (IUNIN,6664) (VOLCOR(I),I=1,idum(0))
          WRITE(IUNOT,6664) (VOLCOR(I),I=1,idum(0))
          deallocate(VOLCOR)
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        ELSE
  241     READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
          IPOS0=INDEX(ZEILE,'*')
          IPOS1=INDEX(ZEILE,'T')
          IPOS2=INDEX(ZEILE,'F')
          IF (IPOS0.EQ.0.AND.IPOS1.EQ.0.AND.IPOS2.EQ.0) GOTO 241
        ENDIF
        IREAD=1
#endif
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
  250   IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        ENDIF
        IF (ZEILE(1:1) .EQ. '*') THEN
          IREAD=0
          GOTO 250
        ENDIF
        READ (ZEILE,6665) ldum(0)
        IREAD=0
#ifndef NEW_SIZE
        READ (IUNIN,6666) idum(0) ! NRADD
        WRITE(IUNOT,6666) idum(0)
        IF (idum(0).GT.0) THEN
          allocate(VOLADD(idum(0)))
          READ (IUNIN,6664) (VOLADD(I),I=1,idum(0))
          WRITE(IUNOT,6664) (VOLADD(I),I=1,idum(0))
          deallocate(VOLADD)
        ENDIF
#else
        IF (ldum(0)) THEN
          READ (IUNIN,6666) idum(0) ! NRADD
          WRITE(IUNOT,6666) idum(0)
          allocate(VOLADD(idum(0)))
          READ (IUNIN,6664) (VOLADD(I),I=1,idum(0))
          WRITE(IUNOT,6664) (VOLADD(I),I=1,idum(0))
          deallocate(VOLADD)
        ELSE
C  FIND START OF NEXT INPUT BLOCK: 3A
  252     READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
          IF (ZEILE(1:3) .NE. '***') GOTO 252
          IREAD=1
        ENDIF    
#endif
*** 3A. DATA FOR NON DEFAULT STANDARD SURFACES
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        ENDIF
  310   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 310
        IREAD=1
        READ (ZEILE,6666) NSTSI
        IREAD=0
        allocate(ILIIN(NLIM+NSTS))
        DO 311 ISTS=1,NSTSI
          IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A72)') ZEILE ! TXTSFL(NLJ)
            WRITE(IUNOT,'(A72)') ZEILE
          ELSE
            IREAD=0
          ENDIF
          READ (IUNIN,6666) idum(0:8) ! JDUMMY,IDIMP,INUMP(ISTS,IDIMP),IRPTA1,
                                      ! IRPTE1,IRPTA2,IRPTE2,IRPTA3,IRPTE3
          IRPTA1 = idum(3)
          IRPTE1 = idum(4) 
          IRPTA2 = idum(5)
          IRPTE2 = idum(6) 
C  OLD INPUT VERSION BEGIN
          IF (idum(1).EQ.1.AND.IRPTA1.NE.IRPTE1) THEN
            idum(3) = convert_eirene
     .       (non_commensurate,2,idum(2),IRPTA1,nx,ny,nx1,ny1,1)
            idum(4) = convert_eirene
     .       (non_commensurate,2,idum(2),IRPTE1,nx,ny,nx1,ny1,3) 
          ELSE
C  NEW INPUT VERSION
            idum(3) = convert_eirene
     .       (non_commensurate,1,IRPTA1,IRPTA2,nx,ny,nx1,ny1,0)
            idum(4) = convert_eirene
     .       (non_commensurate,1,IRPTE1,IRPTE2,nx,ny,nx1,ny1,2)
            idum(5) = convert_eirene
     .       (non_commensurate,2,IRPTA1,IRPTA2,nx,ny,nx1,ny1,1)
            idum(6) = convert_eirene
     .       (non_commensurate,2,IRPTE1,IRPTE2,nx,ny,nx1,ny1,3)
          ENDIF
          WRITE(IUNOT,6666) idum(0:8)
          READ (IUNIN,'(A72)') ZEILE ! ILIIN(NLJ),ILSIDE(NLJ),ILSWCH(NLJ),
                                     ! ILEQUI(NLJ),ILTOR(NLJ),ILCOL(NLJ),
                                     ! ILFIT(NLJ),ILCELL(NLJ),ILBOX(NLJ),
                                     ! ILPLG(NLJ)
          READ (ZEILE,6666) ILIIN (ISTS)
          WRITE(IUNOT,'(A72)') ZEILE
  312     READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
          IREAD=1
          IF (ZEILE(1:1).EQ.'*') THEN
            GOTO 314
#ifdef NEW_SIZE
          ELSEIF (ILIIN (ISTS).LE.0) THEN
            GOTO 312
          ELSEIF (ZEILE(1:8).EQ.'SURFMOD_') THEN
            IREAD=0
            GOTO 314
#endif
          ELSE
            IREAD=0
            READ (IUNIN,'(A72)') ZEILE ! ZNML(NLJ),EWALL(NLJ),EWBIN(NLJ),
                                       ! TRANSP(1,NLJ),TRANSP(2,NLJ),FSHEAT(NLJ)
            WRITE(IUNOT,'(A72)') ZEILE
            READ (IUNIN,'(A72)') ZEILE ! RECYCF(1,NLJ),RECYCT(1,NLJ),RECPRM(1,NLJ),
                                       ! EXPPL(1,NLJ),EXPEL(1,NLJ),EXPIL(1,NLJ)
            WRITE(IUNOT,'(A72)') ZEILE
#ifdef NEW_SIZE
            READ (IUNIN,'(A72)') ZEILE
            WRITE(IUNOT,'(A72)') ZEILE
            IREAD=1
C  READ ONE MORE LINE FOR NON-DEFAULT SPUTTER MODEL
            IF (ZEILE(1:1).NE.'*') IREAD=0
#endif
          ENDIF
  314     CONTINUE
  311   CONTINUE
*** 3B. DATA FOR ADDITIONAL SURFACES           
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        ENDIF
  350   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1).EQ.'*') GOTO 350
        IREAD=0
        READ (ZEILE,6666) NLIMI
        IF (NLIMI.GT.0) THEN
  351     READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
          IREAD=1
          IF (ZEILE(1:3).EQ.'CH0') THEN
            IREAD=0
            GOTO 351
          ENDIF
        ENDIF
        allocate(RLB(NLIMI))
        DO 360 I=1,NLIMI
          IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A72)') ZEILE ! TXTSFL(I)
            WRITE(IUNOT,'(A72)') ZEILE
          ENDIF
          IREAD=0
  362     READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
          IF (ZEILE(1:3).EQ.'CH1') THEN
            GOTO 362
          ELSEIF (ZEILE(1:3).EQ.'CH2') THEN
            GOTO 362
          ELSE
#ifdef NEW_SIZE
            READ (ZEILE,6664) RLB(I),rdum(1:3) ! SAREA(I),RLWMN(I),RLWMX(I)
#else
            READ (ZEILE,6664) RLB(I),rdum(1:2) ! RLWMN(I),RLWMX(I)
#endif
            READ (IUNIN,'(A72)') ZEILE ! ILIIN(I),ILSIDE(I),ILSWCH(I),
                                       ! ILEQUI(I),ILTOR(I),ILCOL(I),
                                       ! ILFIT(I),ILCELL(I),ILBOX(I),
                                       ! ILPLG(I)
            READ (ZEILE,6666) ILIIN(I)
            WRITE(IUNOT,'(A72)') ZEILE
C  READ SURFACE COEFFICIENTS
            IF (RLB(I).LT.2.0_R8) THEN
              READ (IUNIN,'(A72)') ZEILE  ! A0LM(I),A1LM(I),A2LM(I),A3LM(I),A4LM(I),A5LM(I)
              WRITE(IUNOT,'(A72)') ZEILE
              READ (IUNIN,'(A72)') ZEILE  ! A6LM(I),A7LM(I),A8LM(I),A9LM(I)
              WRITE(IUNOT,'(A72)') ZEILE
            ELSEIF (RLB(I).LT.3.0_R8) THEN
              READ (IUNIN,'(A72)') ZEILE  ! P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),P2(3,I)
              WRITE(IUNOT,'(A72)') ZEILE
            ELSEIF (RLB(I).LT.5.0_R8) THEN
              READ (IUNIN,'(A72)') ZEILE  ! P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),P2(3,I)
              WRITE(IUNOT,'(A72)') ZEILE
              READ (IUNIN,'(A72)') ZEILE  ! P3(1,I),P3(2,I),P3(3,I),P4(1,I),P4(2,I),P4(3,I)
              WRITE(IUNOT,'(A72)') ZEILE
            ELSEIF (RLB(I).LT.7.0_R8) THEN
              READ (IUNIN,'(A72)') ZEILE   ! P1(1,I),P1(2,I),P1(3,I),P2(1,I),P2(2,I),P2(3,I)
              WRITE(IUNOT,'(A72)') ZEILE
              READ (IUNIN,'(A72)') ZEILE   ! P3(1,I),P3(2,I),P3(3,I),P4(1,I),P4(2,I),P4(3,I)
              WRITE(IUNOT,'(A72)') ZEILE
              READ (IUNIN,'(A72)') ZEILE   ! P5(1,I),P5(2,I),P5(3,I),P6(1,I),P6(2,I),P6(3,I)
              WRITE(IUNOT,'(A72)') ZEILE
            ENDIF
          ENDIF
C  READ BOUNDARY DATA
          IF (RLB(I).GT.0.0_R8.AND.RLB(I).LT.2.0_R8) THEN
            READ (IUNIN,'(A72)') ZEILE  !  XLIMS1(I,1),YLIMS1(I,1),ZLIMS1(I,1),
                                        !  XLIMS2(I,1),YLIMS2(I,1),ZLIMS2(I,1)
            WRITE(IUNOT,'(A72)') ZEILE
          ELSEIF (RLB(I).LE.0.0_R8) THEN
            IH=-RLB(I)
            idum(1)=B2IDEZ(IH,1,2)
            idum(2)=B2IDEZ(IH,2,2)
            DO 363 J=1,idum(1)
              READ (IUNIN,'(A72)') ZEILE  !  ALIMS(I,J),XLIMS(I,J),YLIMS(I,J),ZLIMS(I,J)
              WRITE(IUNOT,'(A72)') ZEILE
  363       CONTINUE
            DO 364 J=1,idum(2)
              READ (IUNIN,'(A72)') ZEILE  !  ALIMS0(I,J),XLIMS1(I,J),YLIMS1(I,J),
                                          !  ZLIMS1(I,J),XLIMS2(I,J),YLIMS2(I,J)
              WRITE(IUNOT,'(A72)') ZEILE
              READ (IUNIN,'(A72)') ZEILE  !  ZLIMS2(I,J),XLIMS3(I,J),YLIMS3(I,J),ZLIMS3(I,J)
              WRITE(IUNOT,'(A72)') ZEILE
  364       CONTINUE
          ENDIF
C READ LOCAL SURFACE INTERACTION MODEL FOR SURFACE NO. I
  367     READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
          IREAD=1
          IF (ZEILE(1:1).EQ.'*') THEN
            GOTO 369
          ELSEIF (ZEILE(1:9).EQ.'TRANSFORM') THEN
            GOTO 368
          ELSEIF (ILIIN(I).LE.0) THEN
            GOTO 367
          ELSEIF (ZEILE(1:8).EQ.'SURFMOD_') THEN
            IREAD=0
            GOTO 368
          ELSE
            IREAD=0
            READ (IUNIN,'(A72)') ZEILE ! ZNML(I),EWALL(I),EWBIN(I),
                                       ! TRANSP(1,I),TRANSP(2,I),FSHEAT(I)
            WRITE(IUNOT,'(A72)') ZEILE
            READ (IUNIN,'(A72)') ZEILE ! RECYCF(1,I),RECYCT(1,I),RECPRM(1,I),
                                       ! EXPPL(1,I),EXPEL(1,I),EXPIL(1,I)
            WRITE(IUNOT,'(A72)') ZEILE
            READ (IUNIN,'(A72)') ZEILE
            WRITE(IUNOT,'(A72)') ZEILE
            IREAD=1
            IF (ZEILE(1:1).NE.'*'.AND.ZEILE(1:9).NE.'TRANSFORM') IREAD=0
          ENDIF
  368     IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A72)') ZEILE
            WRITE(IUNOT,'(A72)') ZEILE
          ENDIF
          IREAD=1
          IF (ZEILE(1:9).NE.'TRANSFORM') GOTO 369
          IREAD=0
          READ (IUNIN,6666) idum(0:1) ! ITINI,ITEND
          WRITE(IUNOT,6666) idum(0:1)
          READ (IUNIN,'(A72)') ZEILE  ! XLCOR,YLCOR,ZLCOR
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(A72)') ZEILE  ! XLREF,YLREF,ZLREF
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(A72)') ZEILE  ! XLROT,YLROT,ZLROT,ALROT
          WRITE(IUNOT,'(A72)') ZEILE
          GOTO 368
  369     CONTINUE
  360   CONTINUE
        deallocate(RLB)
        deallocate(ILIIN)
C  READ DATA FOR SPECIES SPECIFICATION AND ATOMIC PHYSICS MODULE
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        ENDIF
        IREAD=0
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        READ (IUNIN,*) idum(0)  ! NREACI
        WRITE(IUNOT,6666) idum(0)
  411   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1).EQ.'*') GOTO 420
        READ (ZEILE,66661)
     .      I,FILNAM,H123,REAC,CRC,idum(0:1),rdum(0),RMN,RMX
        IF (RMN.GT.0.0_R8) THEN
          READ (IUNIN,66664) idum(0),(rdum(I),I=1,3)
          WRITE(IUNOT,66664) idum(0),(rdum(I),I=1,3)
        ENDIF
        IF (RMX.GT.0.0_R8) THEN
          READ (IUNIN,66664) idum(0),(rdum(I),I=4,6)
          WRITE(IUNOT,66664) idum(0),(rdum(I),I=4,6)
        ENDIF
        GOTO 411
*4A.   NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES
Cxpb We must read up to section 5 first, to identify ion species
Cxpb and bundles, compute new bulk ion species indices and only
Cxpb then rewrite the entire block
  420   CONTINUE
        READ (IUNIN,*) NATMI   ! last read before possible reassignments
        allocate(TEXTA(NATMI))
        allocate(NMASSA(NATMI),NCHARA(NATMI),NPRTA(NATMI),
     .           ISRFA(NATMI),ISRTA(NATMI),ISRSA(NATMI),
     .           NRCA(NATMI))
#ifdef NEW_SIZE
        allocate(NFOLA(NATMI),NGENA(NATMI))
#endif
        DO 421 IATM=1,NATMI
#ifdef NEW_SIZE
          READ (IUNIN,66666) I,TEXTA(IATM),NMASSA(IATM),NCHARA(IATM),
     .                         NPRTA(IATM),idum(3),
     .                         ISRFA(IATM),ISRTA(IATM),ISRSA(IATM),
     .                         NRCA(IATM),NFOLA(IATM),NGENA(IATM)
#else
          READ (IUNIN,66666) I,TEXTA(IATM),NMASSA(IATM),NCHARA(IATM),
     .                         NPRTA(IATM),idum(3),
     .                         ISRFA(IATM),ISRTA(IATM),ISRSA(IATM),
     .                         NRCA(IATM)
#endif
          DO 422 K=1,NRCA(IATM)
            READ (IUNIN,6666) idum(0:6) 
            READ (IUNIN,6664) rdum(0:4) 
422       CONTINUE
421     CONTINUE
        DO 4211 IATM=1,NATMI
          backspace(IUNIN)
          DO 4221 K=1,NRCA(IATM)
            backspace(IUNIN)
            backspace(IUNIN)
 4221     CONTINUE
 4211   CONTINUE
        NREAC=max(1,maxval(NRCA(1:NATMI)))
        allocate(IREACA(NATMI,NREAC),IBULKA(NATMI,NREAC),
     .           ISCD1A(NATMI,NREAC),ISCD2A(NATMI,NREAC),
     .           ISCDEA(NATMI,NREAC))
#ifdef NEW_SIZE
        allocate(IESTMA(NATMI,NREAC),IBGKA(NATMI,NREAC))
#endif
        allocate(EELECA(NATMI,NREAC),EBULKA(NATMI,NREAC),
     .           ESCD1A(NATMI,NREAC),ESCD2A(NATMI,NREAC),
     .           FREACA(NATMI,NREAC))
        DO 4212 IATM=1,NATMI
          READ (IUNIN,'(A72)') ZEILE
          DO 4222 K=1,NRCA(IATM)
#ifdef NEW_SIZE
            READ (IUNIN,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                        ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                        ISCDEA(IATM,K),IESTMA(IATM,K),
     .                        IBGKA(IATM,K)
#else
            READ (IUNIN,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                        ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                        ISCDEA(IATM,K)
#endif
            READ (IUNIN,6664) EELECA(IATM,K),EBULKA(IATM,K),
     .                        ESCD1A(IATM,K),ESCD2A(IATM,K),
     .                        FREACA(IATM,K)
 4222     CONTINUE
 4212   CONTINUE
*4B.   NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES
        READ (IUNIN,'(A72)') TXTMOL
  430   READ (IUNIN,*) NMOLI
        allocate(TEXTM(NMOLI))
        allocate(NMASSM(NMOLI),NCHARM(NMOLI),NPRTM(NMOLI),
     .           ISRFM(NMOLI),ISRTM(NMOLI),ISRSM(NMOLI),
     .           NRCM(NMOLI))
#ifdef NEW_SIZE
        allocate(NFOLM(NMOLI),NGENM(NMOLI))
#endif
        DO 431 IMOL=1,NMOLI
#ifdef NEW_SIZE
          READ (IUNIN,66666) I,TEXTM(IMOL),NMASSM(IMOL),NCHARM(IMOL),
     .                         NPRTM(IMOL),idum(3),
     .                         ISRFM(IMOL),ISRTM(IMOL),ISRSM(IMOL),
     .                         NRCM(IMOL),NFOLM(IMOL),NGENM(IMOL)
#else
          READ (IUNIN,66666) I,TEXTM(IMOL),NMASSM(IMOL),NCHARM(IMOL),
     .                         NPRTM(IMOL),idum(3),
     .                         ISRFM(IMOL),ISRTM(IMOL),ISRSM(IMOL),
     .                         NRCM(IMOL)
#endif
          DO 432 K=1,NRCM(IMOL)
            READ (IUNIN,6666) idum(0:6) ! IREACM(IMOL,K),IBULKM(IMOL,K),
                                        ! ISCD1M(IMOL,K),ISCD2M(IMOL,K),
                                        ! ISCDEM(IMOL,K),IESTMM(IMOL,K),IBGKM(IMOL,K)
            READ (IUNIN,6664) rdum(0:4) ! EELECM(IMOL,K),EBULKM(IMOL,K),
                                        ! ESCD1M(IMOL,K),ESCD2M(IMOL,K),FREACM(IMOL,K)
432       CONTINUE
431     CONTINUE
        DO 4311 IMOL=1,NMOLI
          backspace(IUNIN)
          DO 4321 K=1,NRCM(IMOL)
            backspace(IUNIN)
            backspace(IUNIN)
 4321     CONTINUE
 4311   CONTINUE
        NREAC=max(1,maxval(NRCM(1:NMOLI)))
        allocate(IREACM(NMOLI,NREAC),IBULKM(NMOLI,NREAC),
     .           ISCD1M(NMOLI,NREAC),ISCD2M(NMOLI,NREAC),
     .           ISCDEM(NMOLI,NREAC))
#ifdef NEW_SIZE
        allocate(IESTMM(NMOLI,NREAC),IBGKM(NMOLI,NREAC))
#endif
        allocate(EELECM(NMOLI,NREAC),EBULKM(NMOLI,NREAC),
     .           ESCD1M(NMOLI,NREAC),ESCD2M(NMOLI,NREAC),
     .           FREACM(NMOLI,NREAC))
        DO 4312 IMOL=1,NMOLI
          READ (IUNIN,'(A72)') ZEILE
          DO 4322 K=1,NRCM(IMOL)
#ifdef NEW_SIZE
            READ (IUNIN,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                        ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                        ISCDEM(IMOL,K),IESTMM(IMOL,K),
     .                        IBGKM(IMOL,K)
#else
            READ (IUNIN,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                        ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                        ISCDEM(IMOL,K)
#endif
            READ (IUNIN,6664) EELECM(IMOL,K),EBULKM(IMOL,K),
     .                        ESCD1M(IMOL,K),ESCD2M(IMOL,K),
     .                        FREACM(IMOL,K)
 4322     CONTINUE
 4312   CONTINUE
*4C.   TEST IONS SPECIES CARDS, NIONI SPECIES 
        READ (IUNIN,'(A72)') TXTION
  440   READ (IUNIN,*) NIONI
        allocate(TEXTI(NIONI))
        allocate(NMASSI(NIONI),NCHARI(NIONI),
     .           NPRTI(NIONI),NCHRGI(NIONI),
     .           ISRFI(NIONI),ISRTI(NIONI),ISRSI(NIONI),
     .           NRCI(NIONI))
#ifdef NEW_SIZE
        allocate(NFOLI(NIONI),NGENI(NIONI))
#endif
        DO 441 IION=1,NIONI
#ifdef NEW_SIZE
          READ (IUNIN,66666) I,TEXTI(IION),NMASSI(IION),NCHARI(IION),
     .                         NPRTI(IION),NCHRGI(IION),
     .                         ISRFI(IION),ISRTI(IION),ISRSI(IION),
     .                         NRCI(IION),NFOLI(IION),NGENI(IION)
#else
          READ (IUNIN,66666) I,TEXTI(IION),NMASSI(IION),NCHARI(IION),
     .                         NPRTI(IION),NCHRGI(IION),
     .                         ISRFI(IION),ISRTI(IION),ISRSI(IION),
     .                         NRCI(IION)
#endif
          DO 442 K=1,NRCI(IION)
            READ (IUNIN,6666) idum(0:6) ! IREACI(IION,K),IBULKI(IION,K),
                                        ! ISCD1I(IION,K),ISCD2I(IION,K),
                                        ! ISCDEI(IION,K),IESTMI(IION,K),IBGKI(IION,K)
            READ (IUNIN,6664) rdum(0:4) ! EELECI(IION,K),EBULKI(IION,K),
                                        ! ESCD1I(IION,K),ESCD2I(IION,K),FREACI(IION,K)
442       CONTINUE
441     CONTINUE
        DO 4411 IION=1,NIONI
          backspace(IUNIN)
          DO 4421 K=1,NRCI(IION)
            backspace(IUNIN)
            backspace(IUNIN)
 4421     CONTINUE
 4411   CONTINUE
        NREAC=max(1,maxval(NRCI(1:NIONI)))
        allocate(IREACI(NIONI,NREAC),IBULKI(NIONI,NREAC),
     .           ISCD1I(NIONI,NREAC),ISCD2I(NIONI,NREAC),
     .           ISCDEI(NIONI,NREAC))
#ifdef NEW_SIZE
        allocate(IESTMI(NIONI,NREAC),IBGKI(NIONI,NREAC))
#endif
        allocate(EELECI(NIONI,NREAC),EBULKI(NIONI,NREAC),
     .           ESCD1I(NIONI,NREAC),ESCD2I(NIONI,NREAC),
     .           FREACI(NIONI,NREAC))
        DO 4412 IION=1,NIONI
          READ (IUNIN,'(A72)') ZEILE
          DO 4422 K=1,NRCI(IION)
#ifdef NEW_SIZE
            READ (IUNIN,6666) IREACI(IION,K),IBULKI(IION,K),
     .                        ISCD1I(IION,K),ISCD2I(IION,K),
     .                        ISCDEI(IION,K),IESTMI(IION,K),
     .                        IBGKI(IION,K)
#else
            READ (IUNIN,6666) IREACI(IION,K),IBULKI(IION,K),
     .                        ISCD1I(IION,K),ISCD2I(IION,K),
     .                        ISCDEI(IION,K)
#endif
            READ (IUNIN,6664) EELECI(IION,K),EBULKI(IION,K),
     .                        ESCD1I(IION,K),ESCD2I(IION,K),
     .                        FREACI(IION,K)
 4422     CONTINUE
 4412   CONTINUE
*** 5. DATA FOR PLASMA BACKGROUND             
        READ (IUNIN,'(A72)') TITLE5
*5A.   BULK ION SPECIES CARDS, NPLSI SPECIES
        READ (IUNIN,'(A72)') TXTPLS
        IREAD = 0
  510   IREAD = IREAD + 1
        call xertst(IREAD.le.100, 'increase size of ZEILEN in b2ytdr!')
        READ (IUNIN,'(A72)') ZEILEN(IREAD)
        IF (ZEILEN(IREAD)(1:1) .EQ. '*') GOTO 510
        READ(ZEILEN(IREAD),6666) NPLSI
        allocate(TEXTP(NPLSI))
        allocate(NMASSP(NPLSI))
        allocate(NCHARP(NPLSI))
        allocate(NPRTP(NPLSI))
        allocate(NCHRGP(NPLSI))
        allocate(ISRFP(NPLSI))
        allocate(ISRTP(NPLSI))
        allocate(ISRSP(NPLSI))
        allocate(NRCP(NPLSI))
        DO 511 IPLS=1,NPLSI
          READ (IUNIN,66666) I,TEXTP(IPLS),NMASSP(IPLS),NCHARP(IPLS),
     .                                     NPRTP(IPLS),NCHRGP(IPLS),
     .                                     ISRFP(IPLS),ISRTP(IPLS),
     .                                     ISRSP(IPLS),NRCP(IPLS)
          DO 512 K=1,NRCP(IPLS)
            READ (IUNIN,6666) idum(0:4) 
            READ (IUNIN,6664) rdum(0:4) 
512       CONTINUE
511     CONTINUE
        DO 5111 IPLS=1,NPLSI
          backspace(IUNIN)
          DO 5121 K=1,NRCP(IPLS)
            backspace(IUNIN)
            backspace(IUNIN)
 5121     CONTINUE
 5111   CONTINUE
        NREAC=max(1,maxval(NRCP(1:NPLSI)))
        allocate(IREACP(NPLSI,NREAC),IBULKP(NPLSI,NREAC),
     .           ISCD1P(NPLSI,NREAC),ISCD2P(NPLSI,NREAC),
     .           ISCDEP(NPLSI,NREAC))
        allocate(EELECP(NPLSI,NREAC),EBULKP(NPLSI,NREAC),
     .           ESCD1P(NPLSI,NREAC),ESCD2P(NPLSI,NREAC),
     .           FREACP(NPLSI,NREAC))
        DO 5112 IPLS=1,NPLSI
          READ (IUNIN,'(A72)') ZEILE
          DO 5122 K=1,NRCP(IPLS)
            READ (IUNIN,6666) IREACP(IPLS,K),IBULKP(IPLS,K),
     .                        ISCD1P(IPLS,K),ISCD2P(IPLS,K),
     .                        ISCDEP(IPLS,K)
            READ (IUNIN,6664) EELECP(IPLS,K),EBULKP(IPLS,K),
     .                        ESCD1P(IPLS,K),ESCD2P(IPLS,K),
     .                        FREACP(IPLS,K)
 5122     CONTINUE
 5112   CONTINUE
cxpb Now we have read the species data and can identify the bundles
        iss = 0
        DO IATM=1,NATMI
          found = .false.
          is = 0
          do while (.not.found .and. is.lt.ns1)
            found = nint(am1(is)).eq.NMASSA(IATM) .and.
     .       nint(zn1(is)).eq.NCHARA(IATM) .and.
     .       is_neutral1(is)
            if (found) then
              iss = iss + 1
              eb2spcr1(iss) = is
              iatm1(iatm) = iss
              in_bundle = .true.
              do while (in_bundle)
                b2espcr1(is) = iss
                is = is + 1
                in_bundle = is.lt.ns1
                if (in_bundle) in_bundle = in_bundle
     &            .and. nint(zn1(is)).eq.NCHARA(IATM)
     &            .and. nint(am1(is)).eq.NMASSA(IATM)
              end do
            else
              is = is + 1
            endif
          enddo
          if (.not.found) then
            iatm1(iatm) = 0
            write(*,*) 'Atom species corresponding to ',
     .       trim(TEXTA(IATM)),' not found! Removing...'
cxpb          call xertst (found, 'Atom species correspondent not found!')
          end if
        END DO
        NATMI1 = iss
        DO 5113 IPLS=1,NPLSI
          found = .false.
          is = 0
          do while (.not.found .and. is.lt.ns1)
            found = nint(am1(is)).eq.NMASSP(IPLS) .and.
     .       nint(zn1(is)).eq.NCHARP(IPLS) .and.
     .       nint(zamin1(is)).le.NCHRGP(IPLS) .and. 
     .       nint(zamax1(is)).ge.NCHRGP(IPLS)
            if (found) then
              eb2ibcr(ipls) = is
            else
              is = is + 1
            endif
          enddo 
          if (.not.found) then
            eb2ibcr(ipls) = -1
            write(*,*) 'Ion species corresponding to ',
     .       trim(TEXTP(IPLS)),' not found! Removing...' 
cxpb          call xertst (found, 'Ion species correspondent not found!')
          end if
 5113   CONTINUE
        e_bundle_start = 0
        e_bundle_end = 0 
        ipls1 = 0 
        ib = 0
        do is = 0, ns1-1
          ipls = 1
          found = .false.
          if (.not.is_neutral1(is)) then
            ib = ib + 1
            do while (IPLS.LE.NPLSI)
              found = eb2ibcr(ipls).eq.is
              if (found) then
                ipls1(ipls) = ib
                if (e_bundle_start(ib).eq.0) then
                  e_bundle_start(ib) = ipls
                  e_bundle_end(ib) = ipls
                else
                  e_bundle_end(ib) = ipls
                endif
              endif
              ipls = ipls + 1
            end do
          end if
        end do
        NPLSI1 = ib
cxpb We now re-compute input data according to the new bundling for block 4
        DO IATM = 1, NATMI
         IF (IATM1(IATM).ne.0) then
          DO K = 1, NRCA(IATM)
            M = B2IDEZ(IBULKA(IATM,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(IBULKA(IATM,K),2,3)
              ICRM = B2IDEZ(IBULKA(IATM,K),3,3)
              IBULKA(IATM,K) = M + 10*L + 100*iatm1(icrm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(IBULKA(IATM,K),2,3)
              IPLS = B2IDEZ(IBULKA(IATM,K),3,3)
              IBULKA(IATM,K) = M + 10*L + 100*ipls1(ipls)
            END IF
            M = B2IDEZ(ISCD1A(IATM,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(ISCD1A(IATM,K),2,3)
              ICRM = B2IDEZ(ISCD1A(IATM,K),3,3)
              ISCD1A(IATM,K) = M + 10*L + 100*iatm1(icrm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(ISCD1A(IATM,K),2,3)
              IPLS = B2IDEZ(ISCD1A(IATM,K),3,3)
              ISCD1A(IATM,K) = M + 10*L + 100*ipls1(ipls)
            END IF
            M = B2IDEZ(ISCD2A(IATM,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(ISCD2A(IATM,K),2,3)
              ICRM = B2IDEZ(ISCD2A(IATM,K),3,3)
              ISCD2A(IATM,K) = M + 10*L + 100*iatm1(icrm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(ISCD2A(IATM,K),2,3)
              IPLS = B2IDEZ(ISCD2A(IATM,K),3,3)
              ISCD2A(IATM,K) = M + 10*L + 100*ipls1(ipls)
            END IF
#ifdef NEW_SIZE
            M = B2IDEZ(IBGKA(IATM,K),1,3)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(IBGKA(IATM,K),2,3)
              ICRM = B2IDEZ(IBGKA(IATM,K),3,3)
              IBGKA(IATM,K) = M + 10*L + 100*iatm1(icrm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(IBGKA(IATM,K),2,3)
              IPLS = B2IDEZ(IBGKA(IATM,K),3,3)
              IBGKA(IATM,K) = M + 10*L + 100*ipls1(ipls)
            END IF
#endif
          END DO
         END IF
        END DO
        DO IMOL = 1, NMOLI
          DO K = 1, NRCM(IMOL)
            M = B2IDEZ(IBULKM(IMOL,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(IBULKM(IMOL,K),2,3)
              IATM = B2IDEZ(IBULKM(IMOL,K),3,3)
              IBULKM(IMOL,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(IBULKM(IMOL,K),2,3)
              IPLS = B2IDEZ(IBULKM(IMOL,K),3,3)
              IBULKM(IMOL,K) = M + 10*L + 100*ipls1(ipls)
            END IF
            M = B2IDEZ(ISCD1M(IMOL,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(ISCD1M(IMOL,K),2,3)
              IATM = B2IDEZ(ISCD1M(IMOL,K),3,3)
              ISCD1M(IMOL,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(ISCD1M(IMOL,K),2,3)
              IPLS = B2IDEZ(ISCD1M(IMOL,K),3,3)
              ISCD1M(IMOL,K) = M + 10*L + 100*ipls1(ipls)
            END IF
            M = B2IDEZ(ISCD2M(IMOL,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(ISCD2M(IMOL,K),2,3)
              IATM = B2IDEZ(ISCD2M(IMOL,K),3,3)
              ISCD2M(IMOL,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(ISCD2M(IMOL,K),2,3)
              IPLS = B2IDEZ(ISCD2M(IMOL,K),3,3)
              ISCD2M(IMOL,K) = M + 10*L + 100*ipls1(ipls)
            END IF
#ifdef NEW_SIZE
            M = B2IDEZ(IBGKM(IMOL,K),1,3)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(IBGKM(IMOL,K),2,3)
              IATM = B2IDEZ(IBGKM(IMOL,K),3,3)
              IBGKM(IMOL,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(IBGKM(IMOL,K),2,3)
              IPLS = B2IDEZ(IBGKM(IMOL,K),3,3)
              IBGKM(IMOL,K) = M + 10*L + 100*ipls1(ipls)
            END IF
#endif
          END DO
        END DO
        DO IION = 1, NIONI
          DO K = 1, NRCI(IION)
            M = B2IDEZ(IBULKI(IION,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(IBULKI(IION,K),2,3)
              IATM = B2IDEZ(IBULKI(IION,K),3,3)
              IBULKI(IION,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(IBULKI(IION,K),2,3)
              IPLS = B2IDEZ(IBULKI(IION,K),3,3)
              IBULKI(IION,K) = M + 10*L + 100*ipls1(ipls)
            END IF
            M = B2IDEZ(ISCD1I(IION,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(ISCD1I(IION,K),2,3)
              IATM = B2IDEZ(ISCD1I(IION,K),3,3)
              ISCD1I(IION,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(ISCD1I(IION,K),2,3)
              IPLS = B2IDEZ(ISCD1I(IION,K),3,3)
              ISCD1I(IION,K) = M + 10*L + 100*ipls1(ipls)
            END IF
            M = B2IDEZ(ISCD2I(IION,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(ISCD2I(IION,K),2,3)
              IATM = B2IDEZ(ISCD2I(IION,K),3,3)
              ISCD2I(IION,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(ISCD2I(IION,K),2,3)
              IPLS = B2IDEZ(ISCD2I(IION,K),3,3)
              ISCD2I(IION,K) = M + 10*L + 100*ipls1(ipls)
            END IF
#ifdef NEW_SIZE
            M = B2IDEZ(IBGKI(IION,K),1,3)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(IBGKI(IION,K),2,3)
              IATM = B2IDEZ(IBGKI(IION,K),3,3)
              IBGKI(IION,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(IBGKI(IION,K),2,3)
              IPLS = B2IDEZ(IBGKI(IION,K),3,3)
              IBGKI(IION,K) = M + 10*L + 100*ipls1(ipls)
            END IF
#endif
          END DO
        END DO
cxpb We now rewrite the modified block 4
*4A.   NEUTRAL ATOMS SPECIES CARDS, NATMI1 SPECIES
        WRITE(IUNOT,6666) NATMI1
        DO 4213 IATM=1,NATMI
         IF (IATM1(IATM).NE.0) THEN
#ifdef NEW_SIZE
          WRITE(IUNOT,66666) IATM1(IATM),
     .                         TEXTA(IATM),NMASSA(IATM),NCHARA(IATM),
     .                         NPRTA(IATM),0,
     .                         ISRFA(IATM),ISRTA(IATM),ISRSA(IATM),
     .                         NRCA(IATM),NFOLA(IATM),NGENA(IATM)
#else
          WRITE(IUNOT,66666) IATM1(IATM),
     .                         TEXTA(IATM),NMASSA(IATM),NCHARA(IATM),
     .                         NPRTA(IATM),0,
     .                         ISRFA(IATM),ISRTA(IATM),ISRSA(IATM),
     .                         NRCA(IATM)
#endif
          DO 4223 K=1,NRCA(IATM)
#ifdef NEW_SIZE
            WRITE(IUNOT,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                        ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                        ISCDEA(IATM,K),IESTMA(IATM,K),
     .                        IBGKA(IATM,K)
#else
            WRITE(IUNOT,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                        ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                        ISCDEA(IATM,K)
#endif
            WRITE(IUNOT,6664) EELECA(IATM,K),EBULKA(IATM,K),
     .                        ESCD1A(IATM,K),ESCD2A(IATM,K),
     .                        FREACA(IATM,K)
 4223     CONTINUE
         END IF
 4213   CONTINUE
*4B.   NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES
        WRITE(IUNOT,'(A72)') TXTMOL
        WRITE(IUNOT,6666) NMOLI
        DO 4313 IMOL=1,NMOLI
#ifdef NEW_SIZE
          WRITE(IUNOT,66666) IMOL,TEXTM(IMOL),NMASSM(IMOL),NCHARM(IMOL),
     .                         NPRTM(IMOL),0,
     .                         ISRFM(IMOL),ISRTM(IMOL),ISRSM(IMOL),
     .                         NRCM(IMOL),NFOLM(IMOL),NGENM(IMOL)
#else
          WRITE(IUNOT,66666) IMOL,TEXTM(IMOL),NMASSM(IMOL),NCHARM(IMOL),
     .                         NPRTM(IMOL),0,
     .                         ISRFM(IMOL),ISRTM(IMOL),ISRSM(IMOL),
     .                         NRCM(IMOL)
#endif
          DO 4323 K=1,NRCM(IMOL)
#ifdef NEW_SIZE
            WRITE(IUNOT,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                        ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                        ISCDEM(IMOL,K),IESTMM(IMOL,K),
     .                        IBGKM(IMOL,K)
#else
            WRITE(IUNOT,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                        ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                        ISCDEM(IMOL,K)
#endif
            WRITE(IUNOT,6664) EELECM(IMOL,K),EBULKM(IMOL,K),
     .                        ESCD1M(IMOL,K),ESCD2M(IMOL,K),
     .                        FREACM(IMOL,K)
 4323     CONTINUE
 4313   CONTINUE
*4C.   TEST IONS SPECIES CARDS, NIONI SPECIES 
        WRITE(IUNOT,'(A72)') TXTION
        WRITE(IUNOT,6666) NIONI
        DO 4413 IION=1,NIONI
#ifdef NEW_SIZE
          WRITE(IUNOT,66666) IION,TEXTI(IION),NMASSI(IION),NCHARI(IION),
     .                         NPRTI(IION),NCHRGI(IION),
     .                         ISRFI(IION),ISRTI(IION),ISRSI(IION),
     .                         NRCI(IION),NFOLI(IION),NGENI(IION)
#else
          WRITE(IUNOT,66666) IION,TEXTI(IION),NMASSI(IION),NCHARI(IION),
     .                         NPRTI(IION),NCHRGI(IION),
     .                         ISRFI(IION),ISRTI(IION),ISRSI(IION),
     .                         NRCI(IION),0,0
#endif
          DO 4423 K=1,NRCI(IION)
#ifdef NEW_SIZE
            WRITE(IUNOT,6666) IREACI(IION,K),IBULKI(IION,K),
     .                        ISCD1I(IION,K),ISCD2I(IION,K),
     .                        ISCDEI(IION,K),IESTMI(IION,K),
     .                        IBGKI(IION,K)
#else
            WRITE(IUNOT,6666) IREACI(IION,K),IBULKI(IION,K),
     .                        ISCD1I(IION,K),ISCD2I(IION,K),
     .                        ISCDEI(IION,K)
#endif
            WRITE(IUNOT,6664) EELECI(IION,K),EBULKI(IION,K),
     .                        ESCD1I(IION,K),ESCD2I(IION,K),
     .                        FREACI(IION,K)
 4423     CONTINUE
 4413   CONTINUE
*** 5. DATA FOR PLASMA BACKGROUND             
        WRITE(IUNOT,'(A72)') TITLE5
*5A.   BULK ION SPECIES CARDS, NPLSI1 SPECIES
        WRITE(IUNOT,'(A72)') TXTPLS
        DO I = 1, IREAD-1
          WRITE(IUNOT,'(A72)') ZEILEN(IREAD)
        END DO
        WRITE(IUNOT,6666) NPLSI1
cxpb Must now convert the bulk ion data according to the new bundling
        allocate(TEXTP1(NPLSI1))
        allocate(NMASSP1(NPLSI1))
        allocate(NCHARP1(NPLSI1))
        allocate(NPRTP1(NPLSI1))
        allocate(NCHRGP1(NPLSI1))
        allocate(ISRFP1(NPLSI1))
        allocate(ISRTP1(NPLSI1))
        allocate(ISRSP1(NPLSI1))
        allocate(NRCP1(NPLSI1))
        allocate(done(0:NPLSI1))
        done = .false.
        DO IS = 1, NPLSI
         IB = IPLS1(IS)
         IF (IB.NE.0 .and. .not.done(IB)) THEN
          IF (e_bundle_start(IB).eq.e_bundle_end(IB)) then
            TEXTP1(IB) = TEXTP(e_bundle_start(IB))
          ELSE
            TEXTP1(IB) = 
     . TEXTP(e_bundle_start(IB))(1:lnblnk(TEXTP(e_bundle_start(IB))))
     . //'_'//
     . TEXTP(e_bundle_end(IB))(1:lnblnk(TEXTP(e_bundle_end(IB))))
          END IF
          NMASSP1(IB) = NMASSP(e_bundle_end(IB))
          NCHARP1(IB) = NCHARP(e_bundle_end(IB))
          NPRTP1(IB) = 0
          NRCP1(IB) = 0
          do ipls = e_bundle_start(IB), e_bundle_end(IB)
            NPRTP1(IB) = NPRTP1(IB) + NPRTP(ipls)
            NRCP1(IB) = NRCP1(IB) + NRCP(ipls)
          end do
          NPRTP1(IB) = 
     .     NPRTP1(IB) / (e_bundle_end(IB)-e_bundle_start(IB)+1)
          NCHRGP1(IB) = 
     .     (NCHRGP(e_bundle_start(IB))+NCHRGP(e_bundle_end(IB)))/2
          ISRFP1(IB) = ISRFP(e_bundle_end(IB))
          ISRTP1(IB) = ISRTP(e_bundle_end(IB))
          ISRSP1(IB) = ISRSP(e_bundle_end(IB))
          done(IB) = .true.
         ENDIF
        ENDDO
        DO 5114 IB=1,NPLSI1
          WRITE(IUNOT,66666) IB,TEXTP1(IB),NMASSP1(IB),NCHARP1(IB),
     .                                     NPRTP1(IB),NCHRGP1(IB),
     .                                     ISRFP1(IB),ISRTP1(IB),
     .                                     ISRSP1(IB),NRCP1(IB)
          DO 5124 IPLS= e_bundle_start(IB), e_bundle_end(IB)
            DO K=1, NRCP(IPLS)
              WRITE(IUNOT,6666) IREACP(IPLS,K),IBULKP(IPLS,K),
     .                          ISCD1P(IPLS,K),ISCD2P(IPLS,K),
     .                          ISCDEP(IPLS,K) 
              WRITE(IUNOT,6664) EELECP(IPLS,K),EBULKP(IPLS,K),
     .                          ESCD1P(IPLS,K),ESCD2P(IPLS,K),
     .                          FREACP(IPLS,K)
            END DO
 5124     CONTINUE
 5114   CONTINUE
        deallocate(NMASSA,NMASSM,NMASSI,NMASSP,NMASSP1)
        deallocate(NCHARA,NCHARM,NCHARI,NCHARP,NCHARP1)
        deallocate(NPRTA,NPRTM,NPRTI,NPRTP,NPRTP1)
        deallocate(NCHRGI,NCHRGP,NCHRGP1)
        deallocate(ISRFA,ISRFM,ISRFI,ISRFP,ISRFP1)
        deallocate(ISRTA,ISRTM,ISRTI,ISRTP,ISRTP1)
        deallocate(ISRSA,ISRSM,ISRSI,ISRSP,ISRSP1)
        deallocate(NRCA,NRCM,NRCI,NRCP,NRCP1)
#ifdef NEW_SIZE
        deallocate(NFOLA,NFOLM,NFOLI)
        deallocate(NGENA,NGENM,NGENI)
#endif
        deallocate(TEXTA,TEXTM,TEXTI,TEXTP,TEXTP1)
        deallocate(IREACA,IBULKA,ISCD1A,ISCD2A,ISCDEA)
        deallocate(EELECA,EBULKA,ESCD1A,ESCD2A,FREACA)
        deallocate(IREACM,IBULKM,ISCD1M,ISCD2M,ISCDEM)
        deallocate(EELECM,EBULKM,ESCD1M,ESCD2M,FREACM)
        deallocate(IREACI,IBULKI,ISCD1I,ISCD2I,ISCDEI)
        deallocate(EELECI,EBULKI,ESCD1I,ESCD2I,FREACI)
        deallocate(IREACP,IBULKP,ISCD1P,ISCD2P,ISCDEP)
        deallocate(EELECP,EBULKP,ESCD1P,ESCD2P,FREACP)
#ifdef NEW_SIZE
        deallocate(IESTMA,IBGKA,IESTMM,IBGKM,IESTMI,IBGKI)
#endif
        IREAD = 0
*5B.   PLASMA BACKGROUND DATA                
  520   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 520
        READ (ZEILE,6666) (INDPRO(J),J=1,12)
        allocate(P0(NPLSI1),P1(NPLSI1),P2(NPLSI1),
     .           P3(NPLSI1),P4(NPLSI1),P5(NPLSI1))
        IF (INDPRO(1).LE.5.AND.NPLSI.GT.0) THEN
          READ (IUNIN,'(A72)') ZEILE ! TE0,TE1,TE2,TE3,TE4,TE5
          WRITE(IUNOT,'(A72)') ZEILE
        END IF
        IF (IABS(INDPRO(2)).LE.5.AND.NPLSI.GT.0) THEN
          IF (INDPRO(2).LT.0) THEN
            allocate(TIO(NPLSI),TIl(NPLSI),TI2(NPLSI),
     .               TI3(NPLSI),TI4(NPLSI),TI5(NPLSI))
            READ (IUNIN,6664) (TIO(I),TIl(I),TI2(I),TI3(I),TI4(I),TI5(I)
     .                        ,I=1,NPLSI)
            DO IB = 1, NPLSI1
              IF (IABS(INDPRO(2)).EQ.1 .OR. IABS(INDPRO(2)).EQ.3) THEN
                P0(IB) = 0.0_R8
                P1(IB) = 0.0_R8
                DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
                  P0(IB) = P0(IB) + TIO(IPLS)
                  P1(IB) = P1(IB) + TIl(IPLS)
                END DO
                P0(IB) = P0(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
                P1(IB) = P1(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
                P2(IB) = TI2(e_bundle_start(IB))
                P3(IB) = TI3(e_bundle_start(IB))
                P4(IB) = TI4(e_bundle_start(IB))
                P5(IB) = TI5(e_bundle_start(IB))              
              ELSE IF (IABS(INDPRO(2)).EQ.2) THEN
                P0(IB) = 0.0_R8
                DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
                  P0(IB) = P0(IB) + TIO(IPLS)
                END DO
                P0(IB) = P0(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
                P1(IB) = TIl(e_bundle_start(IB))
                P2(IB) = TI2(e_bundle_start(IB))
                P3(IB) = TI3(e_bundle_start(IB))
                P4(IB) = TI4(e_bundle_start(IB))
                P5(IB) = TI5(e_bundle_start(IB))
              ELSE
                P0(IB) = TIO(e_bundle_start(IB))
                P1(IB) = TIl(e_bundle_start(IB))
                P2(IB) = TI2(e_bundle_start(IB))
                P3(IB) = TI3(e_bundle_start(IB))
                P4(IB) = TI4(e_bundle_start(IB))
                P5(IB) = TI5(e_bundle_start(IB))
              ENDIF
              WRITE(IUNOT,6664) P0(IB),P1(IB),P2(IB),
     .                          P3(IB),P4(IB),P5(IB)
            ENDDO
            deallocate(TIO,TIl,TI2,TI3,TI4,TI5)
          ELSE
            READ (IUNIN,'(A72)') ZEILE ! TI0(1),TI1(1),TI2(1),TI3(1),TI4(1),TI5(1)
            WRITE(IUNOT,'(A72)') ZEILE
          ENDIF
        ENDIF
        IF (INDPRO(3).LE.5) THEN
          allocate(DI0(NPLSI),DI1(NPLSI),DI2(NPLSI),
     .             DI3(NPLSI),DI4(NPLSI),DI5(NPLSI))
          READ (IUNIN,6664) (DI0(I),DI1(I),DI2(I),DI3(I),DI4(I),DI5(I),
     .                     I=1,NPLSI)
          DO IB = 1, NPLSI1
            IF (IABS(INDPRO(3)).EQ.1 .OR. IABS(INDPRO(3)).EQ.3) THEN
              P0(IB) = 0.0_R8
              P1(IB) = 0.0_R8
              DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
                P0(IB) = P0(IB) + DI0(IPLS)
                P1(IB) = P1(IB) + DI1(IPLS)
              END DO
              P2(IB) = DI2(e_bundle_start(IB))
              P3(IB) = DI3(e_bundle_start(IB))
              P4(IB) = DI4(e_bundle_start(IB))
              P5(IB) = DI5(e_bundle_start(IB))              
            ELSE IF (IABS(INDPRO(3)).EQ.2) THEN
              P0(IB) = 0.0_R8
              DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
                P0(IB) = P0(IB) + DI0(IPLS)
              END DO
              P1(IB) = DI1(e_bundle_start(IB))
              P2(IB) = DI2(e_bundle_start(IB))
              P3(IB) = DI3(e_bundle_start(IB))
              P4(IB) = DI4(e_bundle_start(IB))
              P5(IB) = DI5(e_bundle_start(IB))
            ELSE
              P0(IB) = DI0(e_bundle_start(IB))
              P1(IB) = DI1(e_bundle_start(IB))
              P2(IB) = DI2(e_bundle_start(IB))
              P3(IB) = DI3(e_bundle_start(IB))
              P4(IB) = DI4(e_bundle_start(IB))
              P5(IB) = DI5(e_bundle_start(IB))
            ENDIF
            WRITE(IUNOT,6664) P0(IB),P1(IB),P2(IB),
     .                        P3(IB),P4(IB),P5(IB)
          ENDDO
          deallocate(DI0,DI1,DI2,DI3,DI4,DI5)
        END IF
        IF (IABS(INDPRO(4)).LE.5) THEN
          allocate(VX0(NPLSI),VX1(NPLSI),VX2(NPLSI),
     .             VX3(NPLSI),VX4(NPLSI),VX5(NPLSI))
          READ (IUNIN,6664) (VX0(I),VX1(I),VX2(I),VX3(I),VX4(I),VX5(I),
     .                       I=1,NPLSI)
          DO IB = 1, NPLSI1
            IF (IABS(INDPRO(4)).EQ.1 .OR. IABS(INDPRO(4)).EQ.3) THEN
              P0(IB) = 0.0_R8
              P1(IB) = 0.0_R8
              DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
                P0(IB) = P0(IB) + VX0(IPLS)
                P1(IB) = P1(IB) + VX1(IPLS)
              END DO
              P0(IB) = P0(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
              P1(IB) = P1(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
              P2(IB) = VX2(e_bundle_start(IB))
              P3(IB) = VX3(e_bundle_start(IB))
              P4(IB) = VX4(e_bundle_start(IB))
              P5(IB) = VX5(e_bundle_start(IB))              
            ELSE IF (IABS(INDPRO(4)).EQ.2) THEN
              P0(IB) = 0.0_R8
              DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
                P0(IB) = P0(IB) + VX0(IPLS)
              END DO
              P0(IB) = P0(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
              P1(IB) = VX1(e_bundle_start(IB))
              P2(IB) = VX2(e_bundle_start(IB))
              P3(IB) = VX3(e_bundle_start(IB))
              P4(IB) = VX4(e_bundle_start(IB))
              P5(IB) = VX5(e_bundle_start(IB))
            ELSE
              P0(IB) = VX0(e_bundle_start(IB))
              P1(IB) = VX1(e_bundle_start(IB))
              P2(IB) = VX2(e_bundle_start(IB))
              P3(IB) = VX3(e_bundle_start(IB))
              P4(IB) = VX4(e_bundle_start(IB))
              P5(IB) = VX5(e_bundle_start(IB))
            ENDIF
            WRITE(IUNOT,6664) P0(IB),P1(IB),P2(IB),
     .                        P3(IB),P4(IB),P5(IB)
          ENDDO
          deallocate(VX0,VX1,VX2,VX3,VX4,VX5)
          allocate(VY0(NPLSI),VY1(NPLSI),VY2(NPLSI),
     .             VY3(NPLSI),VY4(NPLSI),VY5(NPLSI))
          READ (IUNIN,6664) (VY0(I),VY1(I),VY2(I),VY3(I),VY4(I),VY5(I),
     .                       I=1,NPLSI)
          DO IB = 1, NPLSI1
            IF (IABS(INDPRO(4)).EQ.1 .OR. IABS(INDPRO(4)).EQ.3) THEN
              P0(IB) = 0.0_R8
              P1(IB) = 0.0_R8
              DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
                P0(IB) = P0(IB) + VY0(IPLS)
                P1(IB) = P1(IB) + VY1(IPLS)
              END DO
              P0(IB) = P0(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
              P1(IB) = P1(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
              P2(IB) = VY2(e_bundle_start(IB))
              P3(IB) = VY3(e_bundle_start(IB))
              P4(IB) = VY4(e_bundle_start(IB))
              P5(IB) = VY5(e_bundle_start(IB))              
            ELSE IF (IABS(INDPRO(4)).EQ.2) THEN
              P0(IB) = 0.0_R8
              DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
                P0(IB) = P0(IB) + VY0(IPLS)
              END DO
              P0(IB) = P0(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
              P1(IB) = VY1(e_bundle_start(IB))
              P2(IB) = VY2(e_bundle_start(IB))
              P3(IB) = VY3(e_bundle_start(IB))
              P4(IB) = VY4(e_bundle_start(IB))
              P5(IB) = VY5(e_bundle_start(IB))
            ELSE
              P0(IB) = VY0(e_bundle_start(IB))
              P1(IB) = VY1(e_bundle_start(IB))
              P2(IB) = VY2(e_bundle_start(IB))
              P3(IB) = VY3(e_bundle_start(IB))
              P4(IB) = VY4(e_bundle_start(IB))
              P5(IB) = VY5(e_bundle_start(IB))
            ENDIF
            WRITE(IUNOT,6664) P0(IB),P1(IB),P2(IB),
     .                        P3(IB),P4(IB),P5(IB)
          ENDDO
          deallocate(VY0,VY1,VY2,VY3,VY4,VY5)
          allocate(VZ0(NPLSI),VZ1(NPLSI),VZ2(NPLSI),
     .             VZ3(NPLSI),VZ4(NPLSI),VZ5(NPLSI))
          READ (IUNIN,6664) (VZ0(I),VZ1(I),VZ2(I),VZ3(I),VZ4(I),VZ5(I),
     .                       I=1,NPLSI)
          DO IB = 1, NPLSI1
            IF (IABS(INDPRO(4)).EQ.1 .OR. IABS(INDPRO(4)).EQ.3) THEN
              P0(IB) = 0.0_R8
              P1(IB) = 0.0_R8
              DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
                P0(IB) = P0(IB) + VZ0(IPLS)
                P1(IB) = P1(IB) + VZ1(IPLS)
              END DO
              P0(IB) = P0(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
              P1(IB) = P1(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
              P2(IB) = VZ2(e_bundle_start(IB))
              P3(IB) = VZ3(e_bundle_start(IB))
              P4(IB) = VZ4(e_bundle_start(IB))
              P5(IB) = VZ5(e_bundle_start(IB))              
            ELSE IF (IABS(INDPRO(4)).EQ.2) THEN
              P0(IB) = 0.0_R8
              DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
                P0(IB) = P0(IB) + VZ0(IPLS)
              END DO
              P0(IB) = P0(IB)/(e_bundle_end(IB)-e_bundle_start(IB)+1)
              P1(IB) = VZ1(e_bundle_start(IB))
              P2(IB) = VZ2(e_bundle_start(IB))
              P3(IB) = VZ3(e_bundle_start(IB))
              P4(IB) = VZ4(e_bundle_start(IB))
              P5(IB) = VZ5(e_bundle_start(IB))
            ELSE
              P0(IB) = VZ0(e_bundle_start(IB))
              P1(IB) = VZ1(e_bundle_start(IB))
              P2(IB) = VZ2(e_bundle_start(IB))
              P3(IB) = VZ3(e_bundle_start(IB))
              P4(IB) = VZ4(e_bundle_start(IB))
              P5(IB) = VZ5(e_bundle_start(IB))
            ENDIF
            WRITE(IUNOT,6664) P0(IB),P1(IB),P2(IB),
     .                        P3(IB),P4(IB),P5(IB)
          ENDDO
          deallocate(VZ0,VZ1,VZ2,VZ3,VZ4,VZ5)
        END IF
        deallocate(P0,P1,P2,P3,P4,P5)
#ifdef NEW_SIZE
        IF (INDPRO(5).EQ.5) THEN
          READ (IUNIN,'(A72)') ZEILE ! BX0,BX1,BX2,BX3,BX4,BX5
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(A72)') ZEILE ! BY0,BY1,BY2,BY3,BY4,BY5
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(A72)') ZEILE ! BZ0,BZ1,BZ2,BZ3,BZ4,BZ5
          WRITE(IUNOT,'(A72)') ZEILE
        ENDIF
#endif
*** 6A. GENERAL DATA FOR REFLECTION MODEL
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
  610   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 610
        IREAD=0
#ifdef NEW_SIZE
        READ (ZEILE,6665) ldum(0) ! NLTRIM
#else
        READ (ZEILE,6665) ldum(0:1) ! NLTRIM,NLMARL
#endif
        IF (ldum(0)) THEN
          READ (IUNIN,'(A72)') ZEILE
          IREAD=1
          IF (INDEX(ZEILE,'PATH')+INDEX(ZEILE,'path').EQ.0) THEN
615         IF (INDEX(ZEILE,'ON')+INDEX(ZEILE,'on').EQ.0) THEN
              GOTO 620
            ELSE
              WRITE(IUNOT,'(A72)') ZEILE
              READ (IUNIN,'(A72)') ZEILE
              GOTO 615
            ENDIF
          ELSE
            WRITE(IUNOT,'(A72)') ZEILE
            READ (IUNIN,'(A72)') ZEILE
  625       IF (INDEX(ZEILE,'ON')+INDEX(ZEILE,'on').EQ.0) THEN
              GOTO 620
            ELSE
              WRITE(IUNOT,'(A72)') ZEILE
              READ (IUNIN,'(A72)') ZEILE
              GOTO 625
            ENDIF
          ENDIF
        ENDIF
 620    allocate(DATM(NATM),DATM1(NATM))
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,6664) (DATM(IATM),IATM=1,NATM)
        ELSE IF (IREAD.EQ.1) THEN
          READ (ZEILE,6664) (DATM(IATM),IATM=1,NATM)
          IREAD=0
        ENDIF
        DATM1(1:NATM) = 0.0_R8
        DO IATM = 1, NATMI
          IF (IATM1(IATM).NE.0) THEN
            DATM1(IATM1(IATM)) = DATM(IATM)
          ENDIF
        ENDDO
        WRITE(IUNOT,6664) DATM1(1:NATM)
        deallocate(DATM,DATM1)
        READ (IUNIN,6664) rdum(1:NMOL) ! (DMOL(IMOL),IMOL=1,NMOL)
        WRITE(IUNOT,6664) rdum(1:NMOL)
        READ (IUNIN,6664) rdum(1:NION) ! (DION(IION),IION=1,NION)
        WRITE(IUNOT,6664) rdum(1:NION)
        allocate(DPLS(NPLS),DPLS1(NPLS))
        READ (IUNIN,6664) (DPLS(IPLS),IPLS=1,NPLS)
        DPLS1(1:NPLS) = 0.0_R8
        DO IB = 1, NPLSI1
          DO IPLS = e_bundle_start(ib), e_bundle_end(ib)
            DPLS1(IB) = DPLS1(IB) + DPLS(IPLS)
          END DO
        END DO
        WRITE(IUNOT,6664) (DPLS1(IPLS),IPLS=1,NPLS)
        deallocate(DPLS,DPLS1)
        READ (IUNIN,'(A72)') ZEILE ! ERMIN,ERCUT,RPROB0,RINTEG,EINTEG,AINTEG
        WRITE(IUNOT,'(A72)') ZEILE
#ifdef NEW_SIZE
        DO
          IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A72)') ZEILE
            WRITE(IUNOT,'(A72)') ZEILE
          END IF
          IF (ZEILE(1:3) .EQ. '***') then
            IREAD = 1
            EXIT
          ENDIF
          IF (ZEILE(1:8) /= 'SURFMOD_') CYCLE
          IREAD=0
          READ (IUNIN,6666) idum(0:3) ! REFCUR%JLREF,REFCUR%JLSPT,
                                      ! REFCUR%JSRS,REFCUR%JSRC
          WRITE(IUNOT,6666) idum(0:3)
          READ (IUNIN,'(A72)') ZEILE ! REFCUR%ZNMLR,REFCUR%EWALLR,REFCUR%EWBINR,
                                     ! REFCUR%TRANSPR(1:2),REFCUR%FSHEATR
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(A72)') ZEILE ! REFCUR%RCYCFR,REFCUR%RCYCTR,REFCUR%RCPRMR,
                                     ! REFCUR%EXPPLR,REFCUR%EXPELR,REFCUR%EXPILR
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
          IREAD=1
          IF (ZEILE(1:1).NE.'*'.AND.ZEILE(1:8).NE.'SURFMOD_') IREAD=0
        END DO
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        ENDIF
#else
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:3) .EQ. '***') GOTO 700
** 6C. DATA FOR CHEMICAL SPUTTERING              
        NCHSPPR=DEF_NATM+DEF_NMOL+DEF_NION+DEF_NPLS
        NCHSPPRI=0
        IREAD=0
  640   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:3) .eq. '***') GOTO 650
        IF (ZEILE(1:1) .eq. '*') GOTO 640
        IF (NCHSPPRI.GE.NCHSPPR) GOTO 640
        NCHSPPRI = NCHSPPRI+1
        READ (IUNIN,6664) rdum(0:5) ! ychsptg(nchsppri),echsptg(nchsppri),
                                    ! fchsptg(nchsppri),(achsptg(i,nchsppri),i=1,3)
        WRITE(IUNOT,6664) rdum(0:5)
        GOTO 640
  650   CONTINUE
#endif
*** 7. DATA FOR PRIMARY SOURCES, NSTRAE STRATA
cxpb  Here must read the whole block before rewrite 
cxpb  to see if some strata need to be removed
  700   CONTINUE
  710   READ (IUNIN,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') THEN
          WRITE(IUNOT,'(A72)') ZEILE
          GOTO 710
        END IF
        READ (ZEILE,6666) NSTRAE
        IREAD=0
        allocate (INDSRC(NSTRAE))
        allocate (TXTSOU(NSTRAE))
        allocate (NLAVRP(NSTRAE),NLAVRT(NSTRAE),
     .            NLSYMP(NSTRAE),NLSYMT(NSTRAE))
        allocate (NPTS(NSTRAE),NINITL(NSTRAE),
     .            NEMODS(NSTRAE),NAMODS(NSTRAE))
        allocate (FLUX(NSTRAE))
        allocate (NLATM(NSTRAE),NLMOL(NSTRAE),
     .            NLION(NSTRAE),NLPLS(NSTRAE))
        allocate (NSPEZ(NSTRAE))
        allocate (NLPNT(NSTRAE),NLLNE(NSTRAE),NLSRF(NSTRAE),
     .            NLVOL(NSTRAE), NLCNS(NSTRAE))
        allocate (NSRFSI(NSTRAE))
        allocate (INDIM(DEF_NSRFS,NSTRAE), INSOR(DEF_NSRFS,NSTRAE),
     .            INGRDA(DEF_NSRFS,NSTRAE,3), 
     .            INGRDE(DEF_NSRFS,NSTRAE,3))
        allocate (SORWGT(DEF_NSRFS,NSTRAE), SORLIM(DEF_NSRFS,NSTRAE),
     .            SORIND(DEF_NSRFS,NSTRAE), SOREXP(DEF_NSRFS,NSTRAE),
     .            SORIFL(DEF_NSRFS,NSTRAE))
        allocate (NRSOR(DEF_NSRFS,NSTRAE), NPSOR(DEF_NSRFS,NSTRAE),
     .            NTSOR(DEF_NSRFS,NSTRAE), NBSOR(DEF_NSRFS,NSTRAE),
     .            NASOR(DEF_NSRFS,NSTRAE), NISOR(DEF_NSRFS,NSTRAE),
     .            ISTOR(DEF_NSRFS,NSTRAE))
        allocate (SORAD1(DEF_NSRFS,NSTRAE), SORAD2(DEF_NSRFS,NSTRAE),
     .            SORAD3(DEF_NSRFS,NSTRAE), SORAD4(DEF_NSRFS,NSTRAE),
     .            SORAD5(DEF_NSRFS,NSTRAE), SORAD6(DEF_NSRFS,NSTRAE))
        allocate (SORENI(NSTRAE), SORENE(NSTRAE), 
     .            SORVDX(NSTRAE), SORVDY(NSTRAE), SORVDZ(NSTRAE))
        allocate (SORCOS(NSTRAE), SORMAX(NSTRAE), 
     .            SORCTX(NSTRAE), SORCTY(NSTRAE), SORCTZ(NSTRAE))
#ifdef NEW_SIZE
        allocate (SCALV(NSTRAE))
        allocate (IVLSF(NSTRAE),
     .            ISCLS(NSTRAE),ISCLT(NSTRAE),ISCL1(NSTRAE),
     .            ISCL2(NSTRAE),ISCL3(NSTRAE),ISCLB(NSTRAE),
     .            ISCLA(NSTRAE))
#endif
        savestratum(0) = .true.
        READ (IUNIN,6666) (INDSRC(ISTRA),ISTRA=1,NSTRAE)
        READ (IUNIN,6664) ALLOC
        DO 712 ISTRA=1,NSTRAE
          IF (INDSRC(ISTRA).EQ.6) THEN
            savestratum(ISTRA) = .true.
            GOTO 712
          ENDIF
          IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A72)') TXTSOU(ISTRA)
          ELSEIF (IREAD.EQ.1) THEN
            READ (ZEILE,'(A72)') TXTSOU(ISTRA)
            IREAD=0
          ENDIF
  713     READ (IUNIN,'(A72)') ZEILE
          IREAD=1
          IF (ZEILE(1:1) .EQ. '*') GOTO 713
          READ (ZEILE,6665) 
     .     NLAVRP(ISTRA),NLAVRT(ISTRA),NLSYMP(ISTRA),NLSYMT(ISTRA)
          IREAD=0
          READ (IUNIN,6666) 
     .     NPTS(ISTRA),NINITL(ISTRA),NEMODS(ISTRA),NAMODS(ISTRA)
#ifndef NEW_SIZE
          READ (IUNIN,6664) FLUX(ISTRA)
#else
          READ (IUNIN,66662) FLUX(ISTRA),SCALV(ISTRA),IVLSF(ISTRA),
     .                       ISCLS(ISTRA),ISCLT(ISTRA),ISCL1(ISTRA),
     .                       ISCL2(ISTRA),ISCL3(ISTRA),ISCLB(ISTRA),
     .                       ISCLA(ISTRA)
          IF (IVLSF(ISTRA).eq.1) THEN ! check volumetric tallies
            IF (ISCLT(ISTRA).eq.1  .or. ISCLT(ISTRA).eq.4 .or. 
     &          ISCLT(ISTRA).eq.8  .or. ISCLT(ISTRA).eq.13 .or. 
     &          ISCLT(ISTRA).eq.18 .or. ISCLT(ISTRA).eq.NTALV-8 .or. 
     &          ISCLT(ISTRA).eq.NTALV-5 .or. 
     &          ISCLT(ISTRA).eq.NTALV-2) THEN ! these are NATM volume tallies
              ISCLS(ISTRA) = iatm1(ISCLS(ISTRA))
            ELSE IF (ISCLT(ISTRA).eq.11 .or. ISCLT(ISTRA).eq.16 .or. 
     &          ISCLT(ISTRA).eq.21) THEN ! these are NPLS volume tallies
              ISCLS(ISTRA) = ipls1(ISCLS(ISTRA))
            END IF
          ELSE IF (IVLSF(ISTRA).eq.2) THEN ! check surface tallies
            IF (ISCLT(ISTRA).eq.1  .or. ISCLT(ISTRA).eq.2  .or. 
     &          ISCLT(ISTRA).eq.3  .or. ISCLT(ISTRA).eq.4  .or.  
     &          ISCLT(ISTRA).eq.17 .or. ISCLT(ISTRA).eq.18 .or. 
     &          ISCLT(ISTRA).eq.19 .or. ISCLT(ISTRA).eq.20 .or.
     &          ISCLT(ISTRA).eq.33) THEN ! these are NATM surface tallies
              ISCLS(ISTRA) = iatm1(ISCLS(ISTRA))
            ELSE IF (ISCLT(ISTRA).eq.36) THEN ! these are NPLS surface tallies
              ISCLS(ISTRA) = ipls1(ISCLS(ISTRA))
            END IF
          END IF
#endif
          READ (IUNIN,6665) 
     .     NLATM(ISTRA),NLMOL(ISTRA),NLION(ISTRA),NLPLS(ISTRA)
          READ (IUNIN,6666) NSPEZ(ISTRA)
          IF (NLPLS(ISTRA)) THEN
            savestratum(ISTRA) = .true.
            IF (NSPEZ(ISTRA).GE.1 .AND. NSPEZ(ISTRA).LE.NPLSI) THEN
              NSPEZ(ISTRA) = ipls1(NSPEZ(ISTRA))
              savestratum(ISTRA) = NSPEZ(ISTRA).NE.0
            ELSE IF (NSPEZ(ISTRA).GT.NPLSI) THEN
              NSPEZ(ISTRA) = max(NSPEZ(ISTRA),NPLSI1+1)
            END IF
          ELSE IF (NLATM(ISTRA)) THEN
            savestratum(ISTRA) = .true.
            IF (NSPEZ(ISTRA).GE.1 .AND. NSPEZ(ISTRA).LE.NATMI) THEN
              NSPEZ(ISTRA) = iatm1(NSPEZ(ISTRA))
              savestratum(ISTRA) = NSPEZ(ISTRA).NE.0
            ELSE IF (NSPEZ(ISTRA).GT.NATMI) THEN
              NSPEZ(ISTRA) = max(NSPEZ(ISTRA),NATMI1+1)
            END IF
          ELSE
            savestratum(ISTRA) = .true.
          END IF
          READ (IUNIN,6665) NLPNT(ISTRA),NLLNE(ISTRA),NLSRF(ISTRA),
     .                      NLVOL(ISTRA),NLCNS(ISTRA)
C  SAME FOR POINT-, LINE-, SURFACE- AND VOLUME SOURCES
          READ (IUNIN,6666) NSRFSI(ISTRA)
          DO 715 J=1,NSRFSI(ISTRA)
            READ (IUNIN,6666) idum(1),INDIM(J,ISTRA),INSOR(J,ISTRA),
     .                        INGRDA(J,ISTRA,1),INGRDE(J,ISTRA,1),
     .                        INGRDA(J,ISTRA,2),INGRDE(J,ISTRA,2),
     .                        INGRDA(J,ISTRA,3),INGRDE(J,ISTRA,3)
            READ (IUNIN,6664) SORWGT(J,ISTRA),SORLIM(J,ISTRA),
     .                        SORIND(J,ISTRA),SOREXP(J,ISTRA),
     .                        SORIFL(J,ISTRA)
            READ (IUNIN,6666) NRSOR(J,ISTRA),NPSOR(J,ISTRA),
     .                        NTSOR(J,ISTRA),NBSOR(J,ISTRA),
     .                        NASOR(J,ISTRA),NISOR(J,ISTRA),
     .                        ISTOR(J,ISTRA)
            READ (IUNIN,6664) SORAD1(J,ISTRA),SORAD2(J,ISTRA),
     .                        SORAD3(J,ISTRA),SORAD4(J,ISTRA),
     .                        SORAD5(J,ISTRA),SORAD6(J,ISTRA)
  715     CONTINUE
C  VELOCITY SPACE DISTRIBUTION
          READ (IUNIN,6664) SORENI(ISTRA),SORENE(ISTRA),
     .                      SORVDX(ISTRA),SORVDY(ISTRA),SORVDZ(ISTRA)
          READ (IUNIN,6664) SORCOS(ISTRA),SORMAX(ISTRA),
     .                      SORCTX(ISTRA),SORCTY(ISTRA),SORCTZ(ISTRA)
  712   CONTINUE
cxpb  Start rewrite of block 7
        IS = 0
        do ISTRA=1,NSTRAE
          if (savestratum(istra)) then
            is = is + 1
            indsrc(is) = indsrc(istra)
          endif
        enddo
        NSTRA1 = is
        savestratum(NSTRAE+1:NSTRA) = .true.
        WRITE(IUNOT,6666) NSTRA1
        WRITE(IUNOT,6666) (INDSRC(ISTRA),ISTRA=1,NSTRA1)
        WRITE(IUNOT,6664) ALLOC
        DO 7121 ISTRA=1,NSTRAE
          IF (INDSRC(ISTRA).EQ.6
     &     .or. .not.savestratum(ISTRA)) GOTO 7121
          WRITE(IUNOT,'(A72)') TXTSOU(ISTRA)
          WRITE(IUNOT,6665) 
     .     NLAVRP(ISTRA),NLAVRT(ISTRA),NLSYMP(ISTRA),NLSYMT(ISTRA)
          WRITE(IUNOT,6666)
     .     NPTS(ISTRA),NINITL(ISTRA),NEMODS(ISTRA),NAMODS(ISTRA)
#ifndef NEW_SIZE
          WRITE(IUNOT,6664) FLUX(ISTRA)
#else
          WRITE(IUNOT,66662) FLUX(ISTRA),SCALV(ISTRA),IVLSF(ISTRA),
     .                       ISCLS(ISTRA),ISCLT(ISTRA),ISCL1(ISTRA),
     .                       ISCL2(ISTRA),ISCL3(ISTRA),ISCLB(ISTRA),
     .                       ISCLA(ISTRA)
#endif
          WRITE(IUNOT,6665) NLATM(ISTRA),NLMOL(ISTRA),
     .                      NLION(ISTRA),NLPLS(ISTRA)
          WRITE(IUNOT,6666) NSPEZ(ISTRA) 
          WRITE(IUNOT,6665) NLPNT(ISTRA),NLLNE(ISTRA),NLSRF(ISTRA),
     .                      NLVOL(ISTRA),NLCNS(ISTRA)
C  SAME FOR POINT-, LINE-, SURFACE- AND VOLUME SOURCES
          WRITE(IUNOT,6666) NSRFSI(ISTRA)
          DO 7151 J=1,NSRFSI(ISTRA)
            idum(2) = INDIM(J,ISTRA)
            idum(3) = INSOR(J,ISTRA)
            idum(8) = INGRDA(J,ISTRA,3)
            idum(9) = INGRDE(J,ISTRA,3)
            IF (INDIM(J,ISTRA).EQ.1.AND.
     &          INGRDA(J,ISTRA,1).NE.INGRDE(J,ISTRA,1)) THEN
C  OLD INPUT VERSION BEGIN
              idum(3) = convert_eirene(non_commensurate,1,
     .         INSOR(J,ISTRA),(INGRDA(J,ISTRA,1)+INGRDE(J,ISTRA,1))/2,
     .          nx,ny,nx1,ny1,4)
              idum(4) = convert_eirene(non_commensurate,2,
     .         INSOR(J,ISTRA),INGRDA(J,ISTRA,1),nx,ny,nx1,ny1,1)
              idum(5) = convert_eirene(non_commensurate,2,
     .         INSOR(J,ISTRA),INGRDE(J,ISTRA,1),nx,ny,nx1,ny1,3)
              idum(6) = INGRDA(J,ISTRA,2)
              idum(7) = INGRDE(J,ISTRA,2)
            ELSE
C  NEW INPUT VERSION
              IF (INDIM(J,ISTRA).EQ.1) THEN
                idum(3) = convert_eirene(non_commensurate,1,
     .           INSOR(J,ISTRA),(INGRDA(J,ISTRA,2)+INGRDE(J,ISTRA,2))/2,
     .           nx,ny,nx1,ny1,4)
              ELSE IF (INDIM(J,ISTRA).EQ.2) THEN
                idum(3) = convert_eirene(non_commensurate,2,
     .           (INGRDA(J,ISTRA,1)+INGRDE(J,ISTRA,1))/2,INSOR(J,ISTRA),
     .           nx,ny,nx1,ny1,4)
              ENDIF
              idum(4) = convert_eirene(non_commensurate,1,
     .         INGRDA(J,ISTRA,1),INGRDA(J,ISTRA,2),nx,ny,nx1,ny1,0)
              idum(5) = convert_eirene(non_commensurate,1,
     .         INGRDE(J,ISTRA,1),INGRDE(J,ISTRA,2),nx,ny,nx1,ny1,2)
              idum(6) = convert_eirene(non_commensurate,2,
     .         INGRDA(J,ISTRA,1),INGRDA(J,ISTRA,2),nx,ny,nx1,ny1,1)
              idum(7) = convert_eirene(non_commensurate,2,
     .         INGRDE(J,ISTRA,1),INGRDE(J,ISTRA,2),nx,ny,nx1,ny1,3)
            ENDIF
            WRITE(IUNOT,6666) J,idum(2:9)
            WRITE(IUNOT,6664) SORWGT(J,ISTRA),SORLIM(J,ISTRA),
     .                        SORIND(J,ISTRA),SOREXP(J,ISTRA),
     .                        SORIFL(J,ISTRA)
            IF (NRSOR(J,ISTRA).GT.0) THEN
              idum(1) = convert_eirene(non_commensurate,1,
     .         NRSOR(J,ISTRA),NPSOR(J,ISTRA),nx,ny,nx1,ny1,4)
            ELSE
              idum(1) = NRSOR(J,ISTRA)
            ENDIF
            IF (NPSOR(J,ISTRA).GT.0) THEN
              idum(2) = convert_eirene(non_commensurate,2,
     .         NRSOR(J,ISTRA),NPSOR(J,ISTRA),nx,ny,nx1,ny1,4)
            ELSE
              idum(2) = NPSOR(J,ISTRA)
            ENDIF
            WRITE(IUNOT,6666) idum(1:2),NTSOR(J,ISTRA),NBSOR(J,ISTRA),
     .                        NASOR(J,ISTRA),NISOR(J,ISTRA),
     .                        ISTOR(J,ISTRA)
            WRITE(IUNOT,6664) SORAD1(J,ISTRA),SORAD2(J,ISTRA),
     .                        SORAD3(J,ISTRA),SORAD4(J,ISTRA),
     .                        SORAD5(J,ISTRA),SORAD6(J,ISTRA)
 7151     CONTINUE
C  VELOCITY SPACE DISTRIBUTION
          WRITE(IUNOT,6664) SORENI(ISTRA),SORENE(ISTRA),
     .                      SORVDX(ISTRA),SORVDY(ISTRA),SORVDZ(ISTRA)
          WRITE(IUNOT,6664) SORCOS(ISTRA),SORMAX(ISTRA),
     .                      SORCTX(ISTRA),SORCTY(ISTRA),SORCTZ(ISTRA)
 7121   CONTINUE
        deallocate (INDSRC)
        deallocate (TXTSOU)
        deallocate (NLAVRP,NLAVRT,NLSYMP,NLSYMT)
        deallocate (NPTS,NINITL,NEMODS,NAMODS)
        deallocate (FLUX)
        deallocate (NLATM,NLMOL,NLION,NLPLS)
        deallocate (NSPEZ)
        deallocate (NLPNT,NLLNE,NLSRF,NLVOL,NLCNS)
        deallocate (NSRFSI)
        deallocate (INDIM,INSOR,INGRDA,INGRDE)
        deallocate (SORWGT,SORLIM,SORIND,SOREXP,SORIFL)
        deallocate (NRSOR,NPSOR,NTSOR,NBSOR,NASOR,NISOR,ISTOR)
        deallocate (SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6)
        deallocate (SORENI,SORENE,SORVDX,SORVDY,SORVDZ)
        deallocate (SORCOS,SORMAX,SORCTX,SORCTY,SORCTZ)
#ifdef NEW_SIZE
        deallocate (SCALV)
        deallocate (IVLSF,ISCLS,ISCLT,ISCL1,ISCL2,ISCL3,ISCLB,ISCLA)
#endif
*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
  810   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 810
        READ (ZEILE,6666) NZADD
        IREAD=0
        DO 811 I=1,NZADD
          IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A72)') ZEILE
            WRITE(IUNOT,'(A72)') ZEILE
          ENDIF
          IREAD=0
          READ (IUNIN,6666) idum(0:1) ! INI,INE
          WRITE(IUNOT,6666) idum(0:1)
  812     READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
          IF (ZEILE(1:3).EQ.'CH3') THEN
            GOTO 812
          ELSEIF (ZEILE(1:1).EQ.'T') THEN
            GOTO 812
          ELSEIF (ZEILE(1:1).EQ.'D') THEN
            GOTO 812
          ELSEIF ((ZEILE(1:1).EQ.'V'.AND.ZEILE(2:2).NE.'L')
     .         .OR.ZEILE(1:1).EQ.'M') THEN
            GOTO 812
          ELSEIF (ZEILE(1:2).EQ.'VL') THEN
            GOTO 812
          ELSEIF (ZEILE(1:1).EQ.'*') THEN
            IREAD=1
            GOTO 811
          ENDIF
  811   CONTINUE
*** 9. DATA FOR STATISTIC AND NONANALOG MODEL
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        ENDIF
        IREAD=0
910     READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') THEN
          WRITE(IUNOT,'(A72)') ZEILE
          GOTO 910
        ELSE
          allocate(NLPRCA(NATM),NLPRCM(NMOL),NLPRCI(NION))
          allocate(NLPRCA1(NATM))
          READ (ZEILE,6665) (NLPRCA(J),J=1,NATM),(NLPRCM(J),J=1,NMOL),
     .                      (NLPRCI(J),J=1,NION) 
          NLPRCA1(1:NATM) = .false.
          DO IATM = 1, NATMI
            IF (IATM1(IATM).NE.0) THEN
              NLPRCA1(IATM1(IATM)) = NLPRCA(IATM)
            ENDIF
          ENDDO
          WRITE(IUNOT,6665) (NLPRCA1(J),J=1,NATM),(NLPRCM(J),J=1,NMOL),
     .                      (NLPRCI(J),J=1,NION) 
          deallocate(NLPRCA,NLPRCM,NLPRCI)
          deallocate(NLPRCA1)
        END IF
#ifndef NEW_SIZE
        NLMAX=MIN0(60,NLIMPS)
        allocate(NSSPL(NLIMPS))
        READ (IUNIN,6665) (ldum(J),J=1,NLMAX) ! (NLPRCF(J),J=1,NLMAX)
        WRITE(IUNOT,6665) (ldum(J),J=1,NLMAX)
        READ (IUNIN,6665) (ldum(J),J=1,NLMAX) ! (NLPRCD(J),J=1,NLMAX)
        WRITE(IUNOT,6665) (ldum(J),J=1,NLMAX)
        READ (IUNIN,6666) idum(0),MAXRAD,(NSSPL(J),J=1,MAXRAD)
        DO J = 1, MAXRAD
          NSSPL(J) = convert_eirene(non_commensurate,
     .     1,NSSPL(J),(nx+1+nncut*3)/2,nx,ny,nx1,ny1,4)
        END DO
        WRITE(IUNOT,6666) idum(0),MAXRAD,(NSSPL(J),J=1,MAXRAD)
        deallocate(NSSPL)
#else
        READ (IUNIN,6666) NPRCSF
        WRITE(IUNOT,6666) NPRCSF
        NPRCSF=MIN0(NLIMPS,NPRCSF)
        IPRCSF=1
  911   CONTINUE
        IF (IPRCSF.LE.NPRCSF) THEN
          READ (IUNIN,6666) (idum(J),J=1,12)
          WRITE(IUNOT,6666) (idum(J),J=1,12)
          IPRCSF=IPRCSF+12
          GOTO 911
        ENDIF
        READ (IUNIN,6666) idum(0),MAXRAD,MAXPOL,MAXTOR,MAXADD
        WRITE(IUNOT,6666) idum(0),MAXRAD,MAXPOL,MAXTOR,MAXADD
        DO I=1,MAXRAD
          READ (IUNIN,66665) idum(0),idum(1),rdum(0) ! ID,NSSPL(IN),PRMSPL(IN)
          idum(1) = convert_eirene(non_commensurate,
     .       1,idum(1),(nx+1+nncut*3)/2,nx,ny,nx1,ny1,4)
          WRITE(IUNOT,66665) idum(0),idum(1),rdum(0)
        ENDDO
        DO I=1,MAXPOL
          READ (IUNIN,66665) idum(0),idum(1),rdum(0)
          idum(1) = convert_eirene(non_commensurate,
     .       2,idum(1),(ny+1)/2,nx,ny,nx1,ny1,4)
          WRITE(IUNOT,66665) idum(0),idum(1),rdum(0)
        ENDDO
        DO I=1,MAXTOR
          READ (IUNIN,66665) idum(0),idum(1),rdum(0)
          WRITE(IUNOT,66665) idum(0),idum(1),rdum(0)
        ENDDO
        DO I=1,MAXADD
          READ (IUNIN,66665) idum(0),idum(1),rdum(0)
          WRITE(IUNOT,66665) idum(0),idum(1),rdum(0)
        ENDDO 
#endif        
C  DATA FOR BIAS SAMPLING
        READ (IUNIN,'(A72)') ZEILE ! WMINV,WMINS,WMINC
        WRITE(IUNOT,'(A72)') ZEILE
        READ (IUNIN,'(A72)') ZEILE ! SPLPAR
        WRITE(IUNOT,'(A72)') ZEILE
C  DATA FOR STANDARD DEVIATION
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
cxpb Read first, transform and rewrite
#ifndef NEW_SIZE
        READ (IUNIN,6666) NSIGVI,NSIGSI,NSIGCI
#else
        READ (IUNIN,6666) NSIGVI,NSIGSI,NSIGCI,NSIGI_BGK,NSIGI_COP
#endif
        allocate (IGH(NSIGVI), IIH(NSIGVI))
        allocate (IGHW(NSIGSI), IIHW(NSIGSI))
        allocate (IGHC(2,NSIGCI), IIHC(2,NSIGCI))
        allocate (savesigv(NSIGVI))
        allocate (savesigs(NSIGSI))
        allocate (savesigc(NSIGCI))
        NSIGV1 = 0
        DO 913 J=1,NSIGVI
          READ (IUNIN,6666) IGH(J),IIH(J)
          IF (IIH(J).eq.1  .or. IIH(J).eq.4  .or. IIH(J).eq.8  .or. 
     &        IIH(J).eq.13 .or. IIH(J).eq.18) THEN ! these are NATM volume tallies
              IGH(J) = iatm1(IGH(J))
#ifdef NEW_SIZE
          ELSE IF (IIH(J).eq.NTALV-8 .or. IIH(J).eq.NTALV-5 .or.
     &        IIH(J).eq.NTALV-2) THEN ! these are NATM volume tallies
              IGH(J) = iatm1(IGH(J))
#endif
          ELSE IF (IIH(J).eq.11 .or. IIH(J).eq.16 .or. IIH(J).eq.21) 
     >      THEN ! these are NPLS volume tallies
              IGH(J) = ipls1(IGH(J))
          END IF
          savesigv(J) = IGH(J).NE.0
          if (savesigv(J)) NSIGV1 = NSIGV1 + 1
  913   CONTINUE
        NSIGS1 = 0
        DO 914 J=1,NSIGSI
          READ (IUNIN,6666) IGHW(J),IIHW(J)
#ifndef NEW_SIZE
          IF (IIHW(J).eq.1  .or. IIHW(J).eq.2  .or. 
     &        IIHW(J).eq.3  .or. IIHW(J).eq.4  .or.  
     &        IIHW(J).eq.13 .or. IIHW(J).eq.14 .or. 
     &        IIHW(J).eq.15 .or. IIHW(J).eq.16 .or.
     &        IIHW(J).eq.25 .or. 
     &        IIHW(J).eq.32 .or. IIHW(J).eq.33) THEN ! these are NATM surface tallies
              IGHW(J) = iatm1(IGHW(J))
          ELSE IF (IIHW(J).eq.28 .or.
     &        IIHW(J).eq.30 .or. IIHW(J).eq.31) THEN ! these are NPLS surface tallies
              IGHW(J) = ipls1(IGHW(J))
          END IF
#else
          IF (IIHW(J).eq.1  .or. IIHW(J).eq.2  .or. 
     &        IIHW(J).eq.3  .or. IIHW(J).eq.4  .or.  
     &        IIHW(J).eq.17 .or. IIHW(J).eq.18 .or. 
     &        IIHW(J).eq.19 .or. IIHW(J).eq.20 .or.
     &        IIHW(J).eq.33) THEN ! these are NATM surface tallies
              IGHW(J) = iatm1(IGHW(J))
          ELSE IF (IIHW(J).eq.36) THEN ! these are NPLS surface tallies
              IGHW(J) = ipls1(IGHW(J))
          END IF
#endif
          savesigs(J) = IGHW(J).NE.0
          if (savesigs(J)) NSIGS1 = NSIGS1 + 1
  914   CONTINUE
        NSIGC1 = 0
        DO 915 J=1,NSIGCI
          READ (IUNIN,6666) IGHC(1,J),IIHC(1,J),IGHC(2,J),IIHC(2,J)
          IF (IIHC(1,J).eq.1  .or. IIHC(1,J).eq.4  .or.  
     &        IIHC(1,J).eq.8  .or. IIHC(1,J).eq.13 .or.  
     &        IIHC(1,J).eq.18) THEN ! these are NATM volume tallies
              IGHC(1,J) = iatm1(IGHC(1,J))
#ifdef NEW_SIZE
          ELSE IF (IIHC(1,J).eq.NTALV-8 .or. IIHC(1,J).eq.NTALV-5 .or.
     &        IIHC(1,J).eq.NTALV-2) THEN ! these are NATM volume tallies
              IGHC(1,J) = iatm1(IGHC(1,J))
#endif
          ELSE IF (IIHC(1,J).eq.11 .or. IIHC(1,J).eq.16 .or.  
     &        IIHC(1,J).eq.21) THEN ! these are NPLS volume tallies
              IGHC(1,J) = ipls1(IGHC(1,J))
          END IF
          IF (IIHC(2,J).eq.1  .or. IIHC(2,J).eq.4  .or.  
     &        IIHC(2,J).eq.8  .or. IIHC(2,J).eq.13 .or.  
     &        IIHC(2,J).eq.18) THEN ! these are NATM volume tallies
              IGHC(2,J) = iatm1(IGHC(2,J))
#ifdef NEW_SIZE
          ELSE IF (IIHC(2,J).eq.NTALV-8 .or. IIHC(2,J).eq.NTALV-5 .or.
     &        IIHC(2,J).eq.NTALV-2) THEN ! these are NATM volume tallies
              IGHC(2,J) = iatm1(IGHC(2,J))
#endif
          ELSE IF (IIHC(2,J).eq.11 .or. IIHC(2,J).eq.16 .or.  
     &        IIHC(2,J).eq.21) THEN ! these are NPLS volume tallies
              IGHC(2,J) = ipls1(IGHC(2,J))
          END IF
          savesigc(J) = IGHC(1,J).NE.0 .AND. IGHC(2,J).NE.0 
          if (savesigc(J)) NSIGC1 = NSIGC1 + 1
  915   CONTINUE
#ifndef NEW_SIZE
        WRITE(IUNOT,6666) NSIGV1,NSIGS1,NSIGC1
#else
        WRITE(IUNOT,6666) NSIGV1,NSIGS1,NSIGC1,NSIGI_BGK,NSIGI_COP
#endif
        DO 9131 J=1,NSIGVI
          if (savesigv(J)) WRITE(IUNOT,6666) IGH(J),IIH(J)
 9131   CONTINUE
        DO 9141 J=1,NSIGSI
          if (savesigs(J)) WRITE(IUNOT,6666) IGHW(J),IIHW(J)
 9141   CONTINUE
        DO 9151 J=1,NSIGCI
          if (savesigc(J)) WRITE(IUNOT,6666) IGHC(1,J),IIHC(1,J),
     .                                       IGHC(2,J),IIHC(2,J)
 9151   CONTINUE
        deallocate (IGH, IIH)
        deallocate (IGHW, IIHW)
        deallocate (IGHC, IIHC)
        deallocate (savesigv)
        deallocate (savesigs)
        deallocate (savesigc)
#ifndef NEW_SIZE
** DATA FOR SOURCE CORRECTION                    
  920   READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3) .EQ. '***') GOTO 1000
        IF (ZEILE(1:1) .EQ. '*') GOTO 920
        K = NATM*NSTRA
        allocate(nlpbls(k),nlpbl1(k))
        IF (K.LE.60) THEN
          READ (ZEILE,6665) (nlpbls(i),i=1,k)
        ELSE
          KMOD = K/60
          DO J = 1,KMOD
            IF (j.eq.1) THEN
              READ (ZEILE,6665) (nlpbls(i),i=(j-1)*60+1,j*60)
            ELSE
              READ (IUNIN,6665) (nlpbls(i),i=(j-1)*60+1,j*60)
            ENDIF
          ENDDO
          READ (IUNIN,6665) (nlpbls(i),i=kmod*60+1,k)
        END IF
        ISTRA = 0
        nlpbl1(1:k) = .false.
        DO j = 1, nstrae
          IF (savestratum(j)) istra = istra+1
          DO i = 1, natmi
            IF (iatm1(i).NE.0 .AND. savestratum(j)) 
     .        nlpbl1(natm*(istra-1)+iatm1(i))=nlpbls(natm*(j-1)+i)
          END DO
        END DO
        WRITE(IUNOT,6665) nlpbl1(1:k)
        deallocate(nlpbls,nlpbl1)
#endif
*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION ESTIMATORS AND ALGEBRAIC EXPRESSIONS
        READ (IUNIN,'(A72)') ZEILE
 1000   WRITE(IUNOT,'(A72)') ZEILE
 1010   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1010
        READ (ZEILE,6666) NADVI,NCLVI,NALVI,NADSI,NALSI
*** 10A. DATA FOR ADDITIONAL TALLIES         
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        DO 1020 J=1,NADVI
 1021     READ (IUNIN,'(A72)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A72)') ZEILE
            GOTO 1021
          ELSE
            READ (ZEILE,6666) idum(0:3) ! IADVE(J),IADVS(J),IADVT(J),IADRC(J)
            IF (idum(2).eq.1  .or. idum(2).eq.4  .or. idum(2).eq.8  .or. 
     &          idum(2).eq.13 .or. idum(2).eq.18) THEN ! these are NATM volume tallies
              idum(1) = iatm1(idum(1))
#ifdef NEW_SIZE
            ELSE IF (idum(2).eq.NTALV-8 .or. idum(2).eq.NTALV-5 .or.
     &        idum(2).eq.NTALV-2) THEN ! these are NATM volume tallies
              idum(1) = iatm1(idum(1))
#endif
            ELSE IF (idum(2).eq.11 .or. idum(2).eq.16 .or. 
     &               idum(2).eq.21) THEN ! these are NPLS volume tallies
              idum(1) = ipls1(idum(1))
            END IF
            WRITE(IUNOT,6666) idum(0:3)
          ENDIF
          READ (IUNIN,'(A72)') ZEILE ! TXTTLA(J)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
          WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 1020   CONTINUE
*** 10B. DATA FOR COLLISION ESTIMATORS
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        DO 1030 J=1,NCLVI
 1031     READ (IUNIN,'(A72)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A72)') ZEILE
            GOTO 1031
          ELSE
            READ (ZEILE,6666) idum(0:3) ! ICLVE(J),ICLVS(J),ICLVT(J),ICLRC(J)
            IF (idum(2).eq.1  .or. idum(2).eq.4  .or. idum(2).eq.8  .or. 
     &          idum(2).eq.13 .or. idum(2).eq.18) THEN ! these are NATM volume tallies
              idum(1) = iatm1(idum(1))
#ifdef NEW_SIZE
          ELSE IF (idum(2).eq.NTALV-8 .or. idum(2).eq.NTALV-5 .or.
     &        idum(2).eq.NTALV-2) THEN ! these are NATM volume tallies
              idum(1) = iatm1(idum(1))
#endif
          ELSE IF (idum(2).eq.11 .or. idum(2).eq.16 .or. idum(2).eq.21)
     >      THEN ! these are NPLS volume tallies
              idum(1) = ipls1(idum(1))
            END IF
            WRITE(IUNOT,6666) idum(0:3)
          ENDIF
          READ (IUNIN,'(A72)') ZEILE ! TXTTLC(J)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
          WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 1030   CONTINUE
*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        DO 1040 M=1,NALVI
 1041     READ (IUNIN,'(A72)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A72)') ZEILE
            GOTO 1041
          ELSE
            READ (ZEILE,'(A72)') CHRTAL
            L = lnblnk(CHRTAL)
            J = INDEX(CHRTAL,',')
            DO WHILE (J.GT.0)
              I = INDEX(CHRTAL,'<')
              K = INDEX(CHRTAL,'>')
              READ(CHRTAL(I+1:J-1),*) idum(1)
              READ(CHRTAL(J+1:K-1),*) idum(2)
!!! must convert <i,j> : species i of volume tally j !!!
              IF (idum(2).eq.1  .or. idum(2).eq.4  .or.  
     &            idum(2).eq.8  .or. idum(2).eq.13 .or. 
     &            idum(2).eq.18) THEN ! these are NATM volume tallies
                idum(1) = iatm1(idum(1))
#ifdef NEW_SIZE
              ELSE IF (idum(2).eq.NTALV-8 .or. idum(2).eq.NTALV-5 .or.
     &            idum(2).eq.NTALV-2) THEN ! these are NATM volume tallies
                idum(1) = iatm1(idum(1))
#endif
              ELSE IF (idum(2).eq.11 .or. idum(2).eq.16 .or. 
     &                 idum(2).eq.21 .or. idum(2).eq.-2 .or.
     &                 idum(2).eq.-4 .or. idum(2).eq.-5 .or.  
     &                 idum(2).eq.-6 .or. idum(2).eq.-7 .or. 
     &                 idum(2).eq.-13) THEN ! these are NPLS volume tallies
                idum(1) = ipls1(idum(1))
              END IF
              IF (idum(1).lt.10) THEN
                WRITE(CHRTAL(I+1:J-1),'(i1)') idum(1)
              ELSE 
                IF (J-1.GT.I+1) then
                  WRITE(CHRTAL(I+1:J-1),'(i2)') idum(1)
                ELSE
                  WRITE(hlp_frm,'(a2,i2,a1)') '(a',L-J+1,')'
                  READ(CHRTAL(J:L),hlp_frm) ZEILE(1:L-J+1)
                  WRITE(CHRTAL(I+1:I+2),'(i2)') idum(1)
                  WRITE(CHRTAL(I+3:I+3+L-J),hlp_frm) ZEILE(1:L-J+1)
                  L = I+3+L-J
                  K = INDEX(CHRTAL(I+3:L),'>') + I+3
                ENDIF
              ENDIF
              J = INDEX(CHRTAL(K+1:L),',')
              IF (J.NE.0) J = J+K
            END DO
            WRITE(IUNOT,'(A72)') CHRTAL
          ENDIF
          READ (IUNIN,'(A72)') ZEILE ! TXTTLR(J)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
          WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 1040   CONTINUE
*** 10D. DATA FOR ADDITIONAL SURFACE TALLIES
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        DO 1050 M=1,NADSI
 1051     READ (IUNIN,'(A72)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A72)') ZEILE
            GOTO 1051
          ELSE
            READ (ZEILE,6666) idum(0:3) ! IADSE(J),IADSS(J),IADST(J),IADSC(J)
#ifndef NEW_SIZE
            IF (idum(2).eq.1  .or. idum(2).eq.2  .or. 
     &          idum(2).eq.3  .or. idum(2).eq.4  .or.  
     &          idum(2).eq.13 .or. idum(2).eq.14 .or. 
     &          idum(2).eq.15 .or. idum(2).eq.16 .or.
     &          idum(2).eq.25 .or. 
     &          idum(2).eq.32 .or. idum(2).eq.33) THEN ! these are NATM surface tallies
              idum(1) = iatm1(idum(1))
            ELSE IF (idum(2).eq.28 .or.
     &          idum(2).eq.30 .or. idum(2).eq.31) THEN ! these are NPLS surface tallies
              idum(1) = ipls1(idum(1))
            END IF
#else
            IF (idum(2).eq.1  .or. idum(2).eq.2  .or. 
     &          idum(2).eq.3  .or. idum(2).eq.4  .or.  
     &          idum(2).eq.17 .or. idum(2).eq.18 .or. 
     &          idum(2).eq.19 .or. idum(2).eq.20 .or.
     &          idum(2).eq.33) THEN ! these are NATM surface tallies
                idum(1) = iatm1(idum(1))
            ELSE IF (idum(2).eq.36) THEN ! these are NPLS surface tallies
                idum(1) = ipls1(idum(1))
            END IF
#endif
            WRITE(IUNOT,6666) idum(0:3)
          ENDIF
          READ (IUNIN,'(A72)') ZEILE ! TXTTLW(J)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
          WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 1050   CONTINUE
*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        DO 1060 M=1,NALSI
 1061     READ (IUNIN,'(A72)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A72)') ZEILE
            GOTO 1061
          ELSE
            READ (ZEILE,'(A72)') CHRTAL
            L = lnblnk(CHRTAL)
            J = INDEX(CHRTAL,',')
            DO WHILE (J.GT.0)
              I = INDEX(CHRTAL,'<')
              K = INDEX(CHRTAL,'>')
              READ(CHRTAL(I+1:J-1),*) idum(1)
              READ(CHRTAL(J+1:K-1),*) idum(2)
!!! must convert <i,j> : species i of surface tally j !!!
#ifndef NEW_SIZE
              IF (idum(2).eq.1  .or. idum(2).eq.2  .or. 
     &            idum(2).eq.3  .or. idum(2).eq.4  .or.  
     &            idum(2).eq.13 .or. idum(2).eq.14 .or. 
     &            idum(2).eq.15 .or. idum(2).eq.16 .or.
     &            idum(2).eq.25 .or. 
     &            idum(2).eq.32 .or. idum(2).eq.33) THEN ! these are NATM surface tallies
                idum(1) = iatm1(idum(1))
              ELSE IF (idum(2).eq.28 .or.
     &            idum(2).eq.30 .or. idum(2).eq.31) THEN ! these are NPLS surface tallies
                idum(1) = ipls1(idum(1))
              END IF
#else
              IF (idum(2).eq.1  .or. idum(2).eq.2  .or. 
     &            idum(2).eq.3  .or. idum(2).eq.4  .or.  
     &            idum(2).eq.17 .or. idum(2).eq.18 .or. 
     &            idum(2).eq.19 .or. idum(2).eq.20 .or.
     &            idum(2).eq.33) THEN ! these are NATM surface tallies
                  idum(1) = iatm1(idum(1))
              ELSE IF (idum(2).eq.36) THEN ! these are NPLS surface tallies
                  idum(1) = ipls1(idum(1))
              END IF
#endif
              IF (idum(1).lt.10) THEN
                WRITE(CHRTAL(I+1:J-1),'(i1)') idum(1)
              ELSE 
                IF (J-1.GT.I+1) then
                  WRITE(CHRTAL(I+1:J-1),'(i2)') idum(1)
                ELSE
                  WRITE(hlp_frm,'(a2,i2,a1)') '(a',L-J+1,')'
                  READ(CHRTAL(J:L),hlp_frm) ZEILE(1:L-J+1)
                  WRITE(CHRTAL(I+1:I+2),'(i2)') idum(1)
                  WRITE(CHRTAL(I+3:I+3+L-J),hlp_frm) ZEILE(1:L-J+1)
                  L = I+3+L-J
                  K = INDEX(CHRTAL(I+3:L),'>') + I+3
                ENDIF
              ENDIF
              J = INDEX(CHRTAL(K+1:L),',')
              IF (J.NE.0) J = J+K
            END DO
            WRITE(IUNOT,'(A72)') CHRTAL
          ENDIF
          READ (IUNIN,'(A72)') ZEILE ! TXTTLW(J)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
          WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 1060   CONTINUE
*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT 
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
 1110   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1110
        allocate(TRCSRC(0:NSTRA),TRCSR1(0:NSTRA))
        READ (IUNIN,6665) (TRCSRC(J),J=0,NSTRA)
        TRCSR1=.false.
        TRCSR1(0)=TRCSRC(0)
        istra = 0
        do J = 1, NSTRA
          if (savestratum(j)) then
            istra = istra+1
            TRCSR1(istra) = TRCSRC(J)
          endif
        enddo
        WRITE(IUNOT,6665) (TRCSR1(J),J=0,NSTRA)
        deallocate(TRCSRC,TRCSR1)
        READ (IUNIN,6666) NVOLPR
        WRITE(IUNOT,6666) NVOLPR
        DO 1120 J=1,NVOLPR
          READ (IUNIN,6666) idum(0:3) ! NTLV,NFLGV,NSPZV1,NSPZV2
          IF (idum(0).eq.1  .or. idum(0).eq.4  .or.  
     &        idum(0).eq.8  .or. idum(0).eq.13 .or. 
     &        idum(0).eq.18) THEN ! these are NATM volume tallies
            idum(4) = iatmbeginrange(idum(2),idum(3),natmi,natmi1,iatm1)
            idum(3) = iatmendrange(idum(2),idum(3),natmi,natmi1,iatm1)
            idum(2) = idum(4)
#ifdef NEW_SIZE
          ELSE IF (idum(0).eq.NTALV-8 .or. idum(0).eq.NTALV-5 .or.
     &        idum(0).eq.NTALV-2) THEN ! these are NATM volume tallies
            idum(4) = iatmbeginrange(idum(2),idum(3),natmi,natmi1,iatm1)
            idum(3) = iatmendrange(idum(2),idum(3),natmi,natmi1,iatm1)
            idum(2) = idum(4)
#endif
          ELSE IF (idum(0).eq.11 .or. idum(0).eq.16 .or. 
     &             idum(0).eq.21 .or. idum(0).eq.-2 .or.
     &             idum(0).eq.-4 .or. idum(0).eq.-5 .or.  
     &             idum(0).eq.-6 .or. idum(0).eq.-7 .or. 
     &             idum(0).eq.-13) THEN ! these are NPLS volume tallies
            idum(4) = iatmbeginrange(idum(2),idum(3),nplsi,nplsi1,ipls1)
            idum(3) = iatmendrange(idum(2),idum(3),nplsi,nplsi1,ipls1)
            idum(2) = idum(4)
          END IF
          WRITE(IUNOT,6666) idum(0:3)
 1120   CONTINUE
        READ (IUNIN,6666) NSURPR
        WRITE(IUNOT,6666) NSURPR
        DO 1130 J=1,NSURPR
#ifndef NEW_SIZE
          READ (IUNIN,6666) idum(0:1) ! NTLS,NFLGS
          WRITE(IUNOT,6666) idum(0:1),0,0
#else
          READ (IUNIN,6666) idum(0:3) ! NTLS,NFLGS,NSPZS1,NSPZS2
          IF (idum(0).eq.1  .or. idum(0).eq.2  .or. 
     &        idum(0).eq.3  .or. idum(0).eq.4  .or.  
     &        idum(0).eq.17 .or. idum(0).eq.18 .or. 
     &        idum(0).eq.19 .or. idum(0).eq.20 .or.
     &        idum(0).eq.33) THEN ! these are NATM surface tallies
            idum(4) = iatmbeginrange(idum(2),idum(3),natmi,natmi1,iatm1)
            idum(3) = iatmendrange(idum(2),idum(3),natmi,natmi1,iatm1)
            idum(2) = idum(4)
          ELSE IF (idum(0).eq.36) THEN ! these are NPLS surface tallies
            idum(4) = iatmbeginrange(idum(2),idum(3),nplsi,nplsi1,ipls1)
            idum(3) = iatmendrange(idum(2),idum(3),nplsi,nplsi1,ipls1)
            idum(2) = idum(4)
          END IF
          WRITE(IUNOT,6666) idum(0:3)
#endif
 1130   CONTINUE
 1131   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1131
#ifndef NEW_SIZE
        READ (IUNIN,6666) NPLINR,NPLOTR,NPLDLR,NPLINP,NPLOTP,NPLDLP,
     .                    NPLINT,NPLOTT,NPLDLT,nsrfcls
#else
        READ (IUNIN,6666) NPLINR,NPLOTR,NPLDLR,NPLINP,NPLOTP,NPLDLP,
     .                    NPLINT,NPLOTT,NPLDLT
#endif        
        idum(0) = convert_eirene
     .   (non_commensurate,1,NPLINR,NPLINP,nx,ny,nx1,ny1,0)
        idum(1) = convert_eirene
     .   (non_commensurate,1,NPLOTR,NPLOTP,nx,ny,nx1,ny1,2)
        idum(2) = max(1, (NPLDLR*(ny1-1))/ny)
        idum(3) = convert_eirene
     .   (non_commensurate,2,NPLINR,NPLINP,nx,ny,nx1,ny1,1)
        idum(4) = convert_eirene
     .   (non_commensurate,2,NPLOTR,NPLOTP,nx,ny,nx1,ny1,3)
        idum(5) = max(1, 
     .   (NPLDLP*(nx1-max(nncut*3,1)))/(nx-max(nncut*3,1)))
#ifndef NEW_SIZE
        WRITE(IUNOT,6666) idum(0:5),NPLINT,NPLOTT,NPLDLT,nsrfcls
        do j=1,abs(nsrfcls) !{
          read (iunin,6666) idum(0) ! msrfcls(j)
          WRITE(IUNOT,6666) idum(0)
          allocate(lsrfcls(idum(0)))
          read (iunin,6666) (lsrfcls(i),i=1,idum(0))
          WRITE(IUNOT,6666) (lsrfcls(i),i=1,idum(0))
          deallocate(lsrfcls)
        end do !}
#else
        WRITE(IUNOT,6666) idum(0:5),NPLINT,NPLOTT,NPLDLT
#endif
C  3D GEOMETRY PLOT
        DO 1140 J=1,5
          READ (IUNIN,6662) ldum(0),TEXTLA,IPLTA
          backspace(IUNIN)
          allocate(IPLAA(IPLTA),IPLEA(IPLTA))
          READ (IUNIN,6662) ldum(0),TEXTLA,IPLTA,
     .                (IPLAA(I),IPLEA(I),I=1,IPLTA)
          WRITE(IUNOT,6662) ldum(0),TEXTLA,IPLTA,
     .                (IPLAA(I),IPLEA(I),I=1,IPLTA)
          deallocate(IPLAA,IPLEA)
 1140   CONTINUE
        DO 1141 J=1,3
          READ (IUNIN,6662) ldum(0),TEXTLA,IPLTA
          backspace(IUNIN)
          allocate(IPLAS(IPLTA),IPLES(IPLTA))
          READ (IUNIN,6662) ldum(0),TEXTLA,IPLTA,
     .                   (IPLAS(I),IPLES(I),I=1,IPLTA)
          IF (J.EQ.1) THEN
            DO I = 1, IPLTA
              IPLAS(I) = convert_eirene
     .         (non_commensurate,1,IPLAS(I),NP2ND/2,nx,ny,nx1,ny1,0)
              IPLES(I) = convert_eirene
     .         (non_commensurate,1,IPLES(I),NP2ND/2,nx,ny,nx1,ny1,2)
            END DO
          ELSE IF (J.EQ.2) THEN
            DO I = 1, IPLTA
              IPLAS(I) = convert_eirene
     .         (non_commensurate,2,NR1ST/2,IPLAS(I),nx,ny,nx1,ny1,1)
              IPLES(I) = convert_eirene
     .         (non_commensurate,2,NR1ST/2,IPLES(I),nx,ny,nx1,ny1,3)
            END DO
          ENDIF
          WRITE(IUNOT,6662) ldum(0),TEXTLA,IPLTA,
     .                (IPLAA(I),IPLEA(I),I=1,IPLTA)
          deallocate(IPLAS,IPLES)
 1141   CONTINUE
        READ (IUNIN,'(A72)') ZEILE ! CH2MX,CH2MY,      CH2X0,CH2Y0,CH2Z0
        WRITE(IUNOT,'(A72)') ZEILE
        READ (IUNIN,'(A72)') ZEILE ! CH3MX,CH3MY,CH3MZ,CH3X0,CH3Y0,CH3Z0
        WRITE(IUNOT,'(A72)') ZEILE
        READ (IUNIN,'(A72)') ZEILE ! ANGLE1,ANGLE2,ANGLE3
        WRITE(IUNOT,'(A72)') ZEILE
        READ (IUNIN,6666) idum(0:10) ! I1TRC,I2TRC,(ISYPLT(J),J=1,8),ILINIE
        WRITE(IUNOT,6666) idum(0:10) 
C  DATA FOR PLOTS OF VOLUME TALLIES
 1151   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1151
        READ (ZEILE,6666) NVOLPL
        IF (NVOLPL.LE.0) GOTO 1175
        DO 1150 J=1,NVOLPL
#ifndef NEW_SIZE
          READ (IUNIN,6666) idum(0),NSP
          WRITE(IUNOT,6666) idum(0),NSP
          READ (IUNIN,6665) ldum(0:3) ! PLTL2D(J),PLTL3D(J),PLTLLG(J),PLTLER(J)
          WRITE(IUNOT,6665) ldum(0:3)
          READ (IUNIN,'(A72)') ZEILE ! TALZMI(J),TALZMA(J)
          WRITE(IUNOT,'(A72)') ZEILE
          IF (ldum(0)) THEN
            READ (IUNIN,6665) ldum(2:3) ! LHIST2(J),LSMOT2(J)
            WRITE(IUNOT,6665) ldum(2:3)
            DO 1160 I=1,NSP
              READ (IUNIN,6666) ISP,NPLIN2,NPLOT2,NPLDL2
              idum(1) = ISP
              IF (idum(0).eq.1  .or. idum(0).eq.4  .or.  
     &            idum(0).eq.8  .or. idum(0).eq.13 .or. 
     &            idum(0).eq.18) THEN ! these are NATM volume tallies
                IF (ISP.GT.0) idum(1) = iatm1(ISP)
              ELSE IF (idum(0).eq.11 .or. idum(0).eq.16 .or. 
     &                 idum(0).eq.21 .or. idum(0).eq.-2 .or.
     &                 idum(0).eq.-4 .or. idum(0).eq.-5 .or.  
     &                 idum(0).eq.-6 .or. idum(0).eq.-7 .or. 
     &                 idum(0).eq.-13) THEN ! these are NPLS volume tallies
                IF (ISP.GT.0) idum(1) = ipls1(ISP)
              END IF
              idum(2) = convert_eirene
     .         (non_commensurate,2,NR1ST/2,NPLIN2,nx,ny,nx1,ny1,1)
              idum(3) = convert_eirene
     .         (non_commensurate,2,NR1ST/2,NPLOT2,nx,ny,nx1,ny1,3)
              idum(4) = max(1, 
     .         (NPLDL2*(nx1-max(nncut*3,1)))/(nx-max(nncut*3,1)))
              WRITE(IUNOT,6666) idum(1:4)
 1160       CONTINUE
          ENDIF
          IF (ldum(1)) THEN
            READ (IUNIN,6665) ldum(2:8) ! LHIST3(J),LCTOUR(J),LSMOT3(J),
                                        ! LRAPS(J),LPRAD(J),LPPOL(J),LPTOR(J)
            WRITE(IUNOT,6665) ldum(2:8)
            DO 1161 I=1,NSP
              READ (IUNIN,6666) ISP, IPROJ, NPLI1, NPLO1, NPLI2, NPLO2
              idum(1) = ISP
              IF (idum(0).eq.1  .or. idum(0).eq.4  .or.  
     &            idum(0).eq.8  .or. idum(0).eq.13 .or. 
     &            idum(0).eq.18) THEN ! these are NATM volume tallies
                IF (ISP.GT.0) idum(1) = iatm1(ISP)
              ELSE IF (idum(0).eq.11 .or. idum(0).eq.16 .or. 
     &                 idum(0).eq.21 .or. idum(0).eq.-2 .or.
     &                 idum(0).eq.-4 .or. idum(0).eq.-5 .or.  
     &                 idum(0).eq.-6 .or. idum(0).eq.-7 .or. 
     &                 idum(0).eq.-13) THEN ! these are NPLS volume tallies
                IF (ISP.GT.0) idum(1) = ipls1(ISP)
              END IF
              IF (ldum(6)) THEN  ! LPRAD
                idum(2) = convert_eirene(non_commensurate,
     .           1,IPROJ,(NPLI1+NPLO1)/2,nx,ny,nx1,ny1,4)
                idum(3) = convert_eirene(non_commensurate,
     .           2,IPROJ,NPLI1,nx,ny,nx1,ny1,1)
                idum(4) = convert_eirene(non_commensurate,
     .           2,IPROJ,NPLO1,nx,ny,nx1,ny1,3)
              ELSE IF (ldum(7)) THEN  ! LPPOL
                idum(2) = convert_eirene(non_commensurate,
     .           2,(NPLI1+NPLO1)/2,IPROJ,nx,ny,nx1,ny1,4)
                idum(3) = convert_eirene(non_commensurate,
     .           1,NPLI1,IPROJ,nx,ny,nx1,ny1,0)
                idum(4) = convert_eirene(non_commensurate,
     .           1,NPLO1,IPROJ,nx,ny,nx1,ny1,2)
              ELSE IF (ldum(8)) THEN  ! LPTOR
                idum(3) = convert_eirene(non_commensurate,
     .           1,NPLI1,NPLI2,nx,ny,nx1,ny1,0)
                idum(4) = convert_eirene(non_commensurate,
     .           1,NPLO1,NPLO2,nx,ny,nx1,ny1,2)
                idum(5) = convert_eirene(non_commensurate,
     .           2,NPLI1,NPLI2,nx,ny,nx1,ny1,1)
                idum(6) = convert_eirene(non_commensurate,
     .           2,NPLO1,NPLO2,nx,ny,nx1,ny1,3)
              ENDIF
              WRITE(IUNOT,6666) idum(1:6)
 1161       CONTINUE
            READ (IUNIN,'(A72)') ZEILE ! TALW1(J),TALW2(J)
            WRITE(IUNOT,'(A72)') ZEILE
          ENDIF
#else
 1152     READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') GOTO 1152
          READ (ZEILE,6666) NSP
          READ (IUNIN,6665) ldum(0:3) ! PLTL2D(J),PLTL3D(J),PLTLLG(J),PLTLER(J)
          WRITE(IUNOT,6665) ldum(0:3)
          READ (IUNIN,'(A72)') ZEILE ! TALZMI(J),TALZMA(J),TALXMI(J),TALXMA(J),
                                     ! TALYMI(J),TALYMA(J)
          WRITE(IUNOT,'(A72)') ZEILE
          IF (ldum(0)) THEN
            READ (IUNIN,6665) ldum(2:3) ! LHIST2(J),LSMOT2(J)
            WRITE(IUNOT,6665) ldum(2:3)
            DO 1160 I=1,NSP
              READ (IUNIN,6666) ISP,NTL,NPLIN2,NPLOT2,NPLDL2
              idum(1) = ISP
              idum(2) = NTL
              IF (idum(2).eq.1  .or. idum(2).eq.4  .or.  
     &            idum(2).eq.8  .or. idum(2).eq.13 .or. 
     &            idum(2).eq.18) THEN ! these are NATM volume tallies
                IF (ISP.GT.0) idum(1) = iatm1(ISP)
              ELSE IF (idum(2).eq.11 .or. idum(2).eq.16 .or. 
     &                 idum(2).eq.21 .or. idum(2).eq.-2 .or.
     &                 idum(2).eq.-4 .or. idum(2).eq.-5 .or.  
     &                 idum(2).eq.-6 .or. idum(2).eq.-7 .or. 
     &                 idum(2).eq.-13) THEN ! these are NPLS volume tallies
                IF (ISP.GT.0) idum(1) = ipls1(ISP)
              END IF
              idum(3) = convert_eirene
     .         (non_commensurate,2,NR1ST/2,NPLIN2,nx,ny,nx1,ny1,1)
              idum(4) = convert_eirene
     .         (non_commensurate,2,NR1ST/2,NPLOT2,nx,ny,nx1,ny1,3)
              idum(5) = max(1, 
     .         (NPLDL2*(nx1-max(nncut*3,1)))/(nx-max(nncut*3,1)))
              WRITE(IUNOT,6666) idum(1:5)
 1160       CONTINUE
          ENDIF
          IF (ldum(1)) THEN
            READ (IUNIN,6665) ldum(2:8) ! LHIST3(J),LCNTR3(J),LSMOT3(J),
                                        ! LRAPS3(J),LVECT3(J),LRPVC3(J)
            WRITE(IUNOT,6665) ldum(2:8)
            READ (IUNIN,6665) ldum(2:4) ! LPRAD3(J),LPPOL3(J),LPTOR3(J)
            WRITE(IUNOT,6665) ldum(2:4)
            DO 1161 I=1,NSP
              READ (IUNIN,6666) ISP,NTL,IPROJ,NPLI1,NPLO1,NPLI2,NPLO2
              idum(1) = ISP
              idum(2) = NTL
              IF (idum(2).eq.1  .or. idum(2).eq.4  .or.  
     &            idum(2).eq.8  .or. idum(2).eq.13 .or. 
     &            idum(2).eq.18) THEN ! these are NATM volume tallies
                IF (ISP.GT.0) idum(1) = iatm1(ISP)
              ELSE IF (idum(2).eq.11 .or. idum(2).eq.16 .or. 
     &                 idum(2).eq.21 .or. idum(2).eq.-2 .or.
     &                 idum(2).eq.-4 .or. idum(2).eq.-5 .or.  
     &                 idum(2).eq.-6 .or. idum(2).eq.-7 .or. 
     &                 idum(2).eq.-13) THEN ! these are NPLS volume tallies
                IF (ISP.GT.0) idum(1) = ipls1(ISP)
              END IF
              IF (ldum(2)) THEN  ! LPRAD
                idum(3) = convert_eirene(non_commensurate,
     .           1,IPROJ,(NPLI1+NPLO1)/2,nx,ny,nx1,ny1,4)
                idum(4) = convert_eirene(non_commensurate,
     .           2,IPROJ,NPLI1,nx,ny,nx1,ny1,1)
                idum(5) = convert_eirene(non_commensurate,
     .           2,IPROJ,NPLO1,nx,ny,nx1,ny1,3)
              ELSE IF (ldum(3)) THEN  ! LPPOL
                idum(3) = convert_eirene(non_commensurate,
     .           2,(NPLI1+NPLO1)/2,IPROJ,nx,ny,nx1,ny1,4)
                idum(4) = convert_eirene(non_commensurate,
     .           1,NPLI1,IPROJ,nx,ny,nx1,ny1,0)
                idum(5) = convert_eirene(non_commensurate,
     .           1,NPLO1,IPROJ,nx,ny,nx1,ny1,2)
              ELSE IF (ldum(4)) THEN  ! LPTOR
                idum(4) = convert_eirene(non_commensurate,
     .           1,NPLI1,NPLI2,nx,ny,nx1,ny1,0)
                idum(5) = convert_eirene(non_commensurate,
     .           1,NPLO1,NPLO2,nx,ny,nx1,ny1,2)
                idum(6) = convert_eirene(non_commensurate,
     .           2,NPLI1,NPLI2,nx,ny,nx1,ny1,1)
                idum(7) = convert_eirene(non_commensurate,
     .           2,NPLO1,NPLO2,nx,ny,nx1,ny1,3)
              ENDIF
              WRITE(IUNOT,6666) idum(1:7)
 1161       CONTINUE
            READ (IUNIN,'(A72)') ZEILE ! TALW1(J),TALW2(J),FCABS1(J),FCABS2(J)
            WRITE(IUNOT,'(A72)') ZEILE
          ENDIF
#endif
 1150   CONTINUE
 1175   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IREAD=1
        IF (ZEILE(1:3) .EQ. '***') GOTO 1205
        GOTO 1175
*** 12. DATA FOR DIAGNOSTIC MODULE
 1205   CONTINUE
 1210   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IREAD=0
        IF (ZEILE(1:1) .EQ. '*') GOTO 1210
#ifndef NEW_SIZE
        READ (ZEILE,6666) NCHORI,idum(0),lntspec,lntsput,idum(1),
     .                    nlstadds,idum(2),ngwties
#else
        READ (ZEILE,6666) NCHORI,idum(0)
#endif
        IF (NCHORI.LE.0) GOTO 1230
        DO 1220 I=1,NCHORI
          READ (IUNIN,'(A72)') ZEILE ! TXTSIG(I)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,6666) idum(0:2) ! NCHTAL(I),NSPSCL(I),NSPNEW(I)
          WRITE(IUNOT,6666) idum(0:2)
#ifndef NEW_SIZE
          READ (IUNIN,6666) idum(0:6) ! NSPSTR(I),NSPSPZ(I),
                                      ! NSPINI(I),NSPEND(I),
                                      ! NSPBLC(I),NSPADD(I),MASLOC
          if (.not.savestratum(idum(0))) idum(0) = 0
          idum(1) = iatm1(idum(1))
          WRITE(IUNOT,6666) idum(0:6)
#else
          READ (IUNIN,6666) idum(0:5) ! NSPSTR(I),NSPSPZ(I),
                                      ! NSPINI(I),NSPEND(I),
                                      ! NSPBLC(I),NSPADD(I)
          if (.not.savestratum(idum(0))) idum(0) = 0
          idum(1) = iatm1(idum(1))
          WRITE(IUNOT,6666) idum(0:5), 0
#endif
          READ (IUNIN,'(A72)') ZEILE ! EMIN1(I),EMAX1(I)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(A72)') ZEILE ! IPIVOT(I),XPILOC,YPILOC,ZPILOC
          WRITE(IUNOT,'(A72)') ZEILE 
          READ (IUNIN,'(A72)') ZEILE ! ICHORD(I),XCHLOC,YCHLOC,ZCHLOC
          WRITE(IUNOT,'(A72)') ZEILE 
 1220   CONTINUE
#ifndef NEW_SIZE
        READ (IUNIN,6665) ldum(0:2) ! PLCHOR,PLSPEC,TRCCHR
        WRITE(IUNOT,6665) ldum(0:2)
 1225   CONTINUE
        IF (lntspec.GT.0) THEN
          READ (IUNIN,'(2e12.0,4i6)') rdum(0:1), idum(0:3) ! emnspc,emxspc,
                                                           ! npspcen,iatspec,npspctt,npspcfi
          WRITE(IUNOT,'(2e12.0,4i6)') rdum(0:1), idum(0:3)
        ENDIF
        DO I=1,lntsput !{
          READ (IUNIN,'(2(i6,2x,a2,2x))') 
     .     idum(0),txtmtmp,idum(1),txtmtra
          WRITE(IUNOT,'(2(i6,2x,a2,2x))') 
     .     idum(0),txtmtmp,idum(1),txtmtra
        END DO !}
        IF (nlstadds.GT.0) THEN !{
          allocate(listadds(nlstadds),listaddg(nlstadds))
          READ (IUNIN,6666) (listadds(i),i=1,nlstadds)
          WRITE(IUNOT,6666) (listadds(i),i=1,nlstadds)
          READ (IUNIN,6666) (listaddg(i),i=1,nlstadds)
          WRITE(IUNOT,6666) (listaddg(i),i=1,nlstadds)
          deallocate(listadds,listaddg)
        END IF !}
        DO j=1,ngwties !{
          READ (IUNIN,6666) rdum(0), mgwties ! kgwties(j),mgwties(j)
          WRITE(IUNOT,6666) rdum(0), mgwties
          IF (mgwties.GT.0) THEN !}{
            allocate(lgwties(mgwties))
            READ (IUNIN,6666) (lgwties(i),i=1,mgwties)
            WRITE(IUNOT,6666) (lgwties(i),i=1,mgwties)
            deallocate(lgwties)
          END IF !}
        END DO !}
 1227   READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IF (ZEILE(1:3).EQ.'CH4') GOTO 1227
        IREAD=1
 1230   CONTINUE
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        END IF
        IREAD=0
        IF (ZEILE(1:1).NE.'*') GOTO 1230
        IREAD=1
#else
        READ (IUNIN,6665) ldum(0:1) ! PLCHOR,PLSPEC
        WRITE(IUNOT,6665) ldum(0:1)
 1230   CONTINUE
C  SKIP READING REST OF THIS BLOCK
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IREAD=0
        IF (ZEILE(1:3).NE.'***') GOTO 1230
        IREAD=1
#endif
*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        END IF
        IREAD=0
        READ (IUNIN,6666) idum(0) ! NPRNLI
        WRITE(IUNOT,6666) idum(0)
        IF (idum(0).LE.0) GOTO 1350
        READ (IUNIN,6666) idum(0:1) ! NTIMV,NPTST
        WRITE(IUNOT,6666) idum(0:1)
        READ (IUNIN,'(A72)') ZEILE ! DTIMV,TIME0
        WRITE(IUNOT,'(A72)') ZEILE
*** 13A. DATA FOR SNAPSHOT TALLIES           
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        READ (IUNIN,6666) NSNVI
        WRITE(IUNOT,6666) NSNVI
        DO 1320 J=1,NSNVI
 1321     READ (IUNIN,'(A72)') ZEILE ! ISNVE(J),ISNVS(J),ISNVT(J),ISNRC(J)
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A72)') ZEILE
            GOTO 1321
          ELSE
            READ (ZEILE,6666) idum(0:3) ! ISNVE(J),ISNVS(J),ISNVT(J),ISNRC(J)
            IF (idum(2).eq.1  .or. idum(2).eq.4  .or. idum(2).eq.8  .or. 
     &          idum(2).eq.13 .or. idum(2).eq.18) THEN ! these are NATM volume tallies
              idum(1) = iatm1(idum(1))
#ifdef NEW_SIZE
            ELSE IF (idum(2).eq.NTALV-8 .or. idum(2).eq.NTALV-5 .or.
     &        idum(2).eq.NTALV-2) THEN ! these are NATM volume tallies
              idum(1) = iatm1(idum(1))
#endif
            ELSE IF (idum(2).eq.11 .or. idum(2).eq.16 .or.  
     &               idum(2).eq.21) THEN ! these are NPLS volume tallies
              idum(1) = ipls1(idum(1))
            END IF
            WRITE(IUNOT,6666) idum(0:3) ! ISNVE(J),ISNVS(J),ISNVT(J),ISNRC(J)
          ENDIF         
          READ (IUNIN,'(A72)') ZEILE ! TXTTAL(J,NTALT)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
          WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 1320   CONTINUE
C  SKIP READING REST OF THIS BLOCK
 1350   CONTINUE
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        IREAD=0
        IF (ZEILE(1:3).NE.'***') GOTO 1350
        IREAD=1
*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        END IF
        IF (NMODE.EQ.0) THEN
#ifdef NEW_SIZE
          READ (IUNIN,6666) NAINI
          WRITE(IUNOT,6666) NAINI
          DO 3020 J=1,NAINI
 3021       READ (IUNIN,'(A72)') ZEILE
            WRITE(IUNOT,'(A72)') ZEILE
            IF (ZEILE(1:1) .EQ. '*') GOTO 3021
            READ (IUNIN,6666) idum(0:1) ! NAINS(J),NAINT(J)
            IF (idum(1).EQ.22 .or. idum(1).EQ.122) THEN ! multiple NPLS arrays
              idum(2) = idum(0)/NPLS
              idum(3) = idum(0) - NPLS*idum(2)
              idum(3) = ipls1(idum(3))
              idum(0) = idum(3) + NPLS*idum(2)
            ELSE IF (idum(1).GE.17) THEN ! Eirene arrays
              idum(0) = ipls1(idum(0))
            ELSE  ! B2 arrays
              idum(0) = is1(idum(0))
            ENDIF
            WRITE(IUNOT,6666) idum(0:1)
            READ (IUNIN,'(A72)') ZEILE ! TXTPLS(J,NTALN)
            WRITE(IUNOT,'(A72)') ZEILE
            READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
            WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 3020     CONTINUE
#endif
        ELSE IF (NMODE.NE.0) THEN
          READ (IUNIN,'(A72)') ZEILE ! LSYMET,LBALAN
          WRITE(IUNOT,'(A72)') ZEILE
#ifndef NEW_SIZE
          READ (IUNIN,'(A72)') ZEILE ! NDUMMY
          WRITE(IUNOT,'(A72)') ZEILE
#endif
          READ (IUNIN,'(3I6)') idum(0:2) ! NFLA,NCUTB,NCUTL
          idum(0) = ns1-nspecies1
          WRITE(IUNOT,'(3I6)') idum(0:2)
          allocate(IFLB(NPLSI),FCTE(NPLSI),BMASS(NPLSI))
          allocate(IFLB1(NPLSI1),FCTE1(NPLSI1),BMASS1(NPLSI1))
#ifndef NEW_SIZE
          allocate(lkindp(NPLSI))
          allocate(lkindp1(NPLSI1))
#endif
          DO I=1,NPLSI          
#ifndef NEW_SIZE
            READ (IUNIN,'(2I6,2E12.4,I6)')
     .       idum(0),IFLB(I),FCTE(I),BMASS(I),LKINDP(I)
            IF (ipls1(i).NE.0) lkindp1(ipls1(i)) = lkindp(i)
#else
            READ (IUNIN,'(2I6,2E12.4)') idum(0),IFLB(I),FCTE(I),BMASS(I)
#endif            
            IF (ipls1(i).NE.0) THEN
              IFLB1(ipls1(I)) = ifl1(eb2ibcr(I))
              FCTE1(ipls1(I)) = FCTE(I)
              BMASS1(ipls1(I)) = BMASS(I)
            END IF
          END DO
          DO I=1,NPLSI1
#ifndef NEW_SIZE
            WRITE(IUNOT,'(1P,2I6,2E12.4,I6)')
     .       I,IFLB1(I),FCTE1(I),BMASS1(I),LKINDP1(I)
#else
            WRITE(IUNOT,'(1P,2I6,2E12.4)') I,IFLB1(I),FCTE1(I),BMASS1(I)
#endif
          END DO
          deallocate(IFLB,FCTE,BMASS)
          deallocate(IFLB1,FCTE1,BMASS1)
#ifndef NEW_SIZE
          deallocate(lkindp)
          deallocate(lkindp1)
#endif
          READ (IUNIN,'(2I6)') idum(0:1) ! NDXA,NDYA
          WRITE(IUNOT,'(2I6)') nx1+max(1,nncut*3)-1, ny1
C  NUMBER OF TARGETS ON B2 SURFACES: NTARGI
          READ (IUNIN,'(I6)') NTARGI
          WRITE(IUNOT,'(I6)') NTARGI
          allocate(NTGPRT(NTARGI))
#ifdef NEW_SIZE
          IF (NTARGI.GT.0) THEN
            READ (IUNIN,6666) (NTGPRT(I),I=1,NTARGI)
            WRITE(IUNOT,6666) (NTGPRT(I),I=1,NTARGI)
          END IF
#else
          NTGPRT(1:NTARGI) = 1
#endif
          DO 3030 I=1,NTARGI
            DO 3033 J=1,NTGPRT(I)
#ifndef NEW_SIZE
              READ (IUNIN,6666) idum(0),NDT,NINCT,NIXY,NTIN,
     .                          NTEN,NIFLG,NPTC
#else
              READ (IUNIN,6666) idum(0),NDT,NINCT,NIXY,NTIN,
     .                          NTEN,NIFLG,NPTC,NSPZI,NSPZE
              idum(8) = is1beginrange(NSPZI,NSPZE,ns,ns1,is1)
              idum(9) = is1endrange(NSPZI,NSPZE,ns,ns1,is1)
#endif
              IF (NIXY.EQ.1) THEN  ! target runs radially
                IF (NINCT.EQ.1) THEN
                  iside = 3
                ELSE IF (NINCT.EQ.-1) THEN
                  iside = 1
                ENDIF
                idum(1) = convert_eirene(non_commensurate,
     .           2,(NTIN+NTEN)/2,NDT,nx,ny,nx1,ny1,iside)
                idum(4) = convert_eirene(non_commensurate,
     .           1,NTIN,NDT,nx,ny,nx1,ny1,0)
                idum(5) = convert_eirene(non_commensurate,
     .           1,NTEN,NDT,nx,ny,nx1,ny1,2)
              ELSE IF (NIXY.EQ.2) THEN  ! target runs poloidally
                IF (NINCT.EQ.1) THEN
                  iside = 2
                ELSE IF (NINCT.EQ.-1) THEN
                  iside = 0
                ENDIF
                idum(1) = convert_eirene(non_commensurate,
     .           1,NDT,(NTIN+NTEN)/2,nx,ny,nx1,ny1,iside)
                idum(4) = convert_eirene(non_commensurate,
     .           2,NDT,NTIN,nx,ny,nx1,ny1,1)
                idum(5) = convert_eirene(non_commensurate,
     .           2,NDT,NTEN,nx,ny,nx1,ny1,3)
              ENDIF
#ifndef NEW_SIZE
              WRITE(IUNOT,6666) idum(0:1),NINCT,NIXY,
     .         idum(4:5),NIFLG,NPTC
#else
              WRITE(IUNOT,6666) idum(0:1),NINCT,NIXY,
     .         idum(4:5),NIFLG,NPTC,idum(8:9)
#endif
 3033       CONTINUE
 3030     CONTINUE
          READ (IUNIN,'(A72)') ZEILE ! CHGP,CHGEE,CHGEI,CHGMOM
          WRITE(IUNOT,'(A72)') ZEILE
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2 INTO EIRENE
C  HERE: B2 VOLUME TALLIES
          READ (IUNIN,'(I6)') NAINB
          WRITE(IUNOT,'(I6)') NAINB
          DO 3040 I=1,NAINB
            READ (IUNIN,'(6I6)') idum(0:2) ! I,NAINS(I),NAINT(I)
            IF (idum(2).EQ.22 .or. idum(2).EQ.122) THEN ! multiple NPLS arrays
              idum(3) = idum(1)/NPLS
              idum(4) = idum(1) - NPLS*idum(3)
              idum(4) = ipls1(idum(4))
              idum(1) = idum(4) + NPLS*idum(3)
            ELSE IF (idum(2).GE.17) THEN ! Eirene arrays
              idum(1) = ipls1(idum(1))
            ELSE  ! B2 arrays
              idum(1) = is1(idum(1))
            ENDIF
            WRITE(IUNOT,'(6I6)') idum(0:2)
            READ (IUNIN,'(A72)') ZEILE ! TXTPLS(I,11)
            WRITE(IUNOT,'(A72)') ZEILE
            READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
            WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 3040     CONTINUE
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO B2
C  HERE: EIRENE SURFACE TALLIES
          READ (IUNIN,'(I6)') NAOTB
          WRITE(IUNOT,'(I6)') NAOTB
          DO 3050 I=1,NAOTB
            READ (IUNIN,'(6I6)') idum(0:2) ! I,NAOTS(IAOT),NAOTT(IAOT)
#ifndef NEW_SIZE
            IF (idum(2).eq.1  .or. idum(2).eq.2  .or. 
     &          idum(2).eq.3  .or. idum(2).eq.4  .or.  
     &          idum(2).eq.13 .or. idum(2).eq.14 .or. 
     &          idum(2).eq.15 .or. idum(2).eq.16 .or.
     &          idum(2).eq.25 .or. 
     &          idum(2).eq.32 .or. idum(2).eq.33) THEN ! these are NATM surface tallies
              idum(1) = iatm1(idum(1))
            ELSE IF (idum(2).eq.28 .or.
     &          idum(2).eq.30 .or. idum(2).eq.31) THEN ! these are NPLS surface tallies
              idum(1) = ipls1(idum(1))
            END IF
#else
            IF (idum(2).eq.1  .or. idum(2).eq.2  .or. 
     &          idum(2).eq.3  .or. idum(2).eq.4  .or.  
     &          idum(2).eq.17 .or. idum(2).eq.18 .or. 
     &          idum(2).eq.19 .or. idum(2).eq.20 .or.
     &          idum(2).eq.33) THEN ! these are NATM surface tallies
                idum(1) = iatm1(idum(1))
            ELSE IF (idum(2).eq.36) THEN ! these are NPLS surface tallies
                idum(1) = ipls1(idum(1))
            END IF
#endif
            WRITE(IUNOT,'(6I6)') idum(0:2)
 3050     CONTINUE
        ENDIF
*** 15. ADDITIONAL GEOMETRY
        READ (IUNIN,'(a80)',end=8686) geometry_comment
        WRITE(IUNOT,'(a80)') geometry_comment
        IF (NPPLG.EQ.1) THEN             ! limiter case
          numpairs=1
        ELSE IF (NPPLG.EQ.3) THEN         ! single null
          numpairs=2
        ELSE IF (NPPLG.EQ.6) THEN         ! double null
          numpairs=4
        END IF
        DO I=1,numpairs*2
          READ (IUNIN,*) idum(0:1) ! onetwo(i),limpos(i)
          WRITE(IUNOT,6666) idum(0:1)
          IF (idum(0).LT.0) THEN
            READ (IUNIN,*) rdum(0:1) ! xpolpos(i),ypolpos(i)
            WRITE(IUNOT,6664) rdum(0:1)
          END IF
        END DO
 8686   continue
        close(IUNOT)
        goto 8688
 8687   continue
        write(*,*) 'Empty input.dat file not converted'
        savestratum(0:nstrai) = .true.
 8688   continue
        close(IUNIN)
      else
        savestratum(0:nstrai) = .true.
      endif 

*  ..rewrite neutrals namelist
*  ..recall it was already read earlier
      inquire(file='../b2.neutrals.parameters',exist=todo)
      if (neutrals_namelist.eq.1.or.use_eirene.ne.0.or.todo) then
        ib1 = 0
        do ib = 1, nstrai
         if (savestratum(ib) .or. crcstra(ib).eq.'T') then
          ib1 = ib1+1
          crcstra1(ib1) = crcstra(ib)
          if (crcstra(ib).eq.'S') then
            rcpos1(ib1) = conformal_convert(non_commensurate,
     .       2,rcstart(ib),rcpos(ib),nx,ny,nx1,ny1,2)
            rcstart1(ib1) = conformal_convert(non_commensurate,
     .       1,rcstart(ib),rcpos(ib),nx,ny,nx1,ny1,1) 
            rcend1(ib1) = conformal_convert(non_commensurate,
     .       1,rcend(ib),rcpos(ib),nx,ny,nx1,ny1,3)
          else if (crcstra(ib).eq.'N') then
            rcpos1(ib1) = conformal_convert(non_commensurate,
     .       2,rcstart(ib),rcpos(ib),nx,ny,nx1,ny1,0)
            rcstart1(ib1) = conformal_convert(non_commensurate,
     .       1,rcstart(ib),rcpos(ib),nx,ny,nx1,ny1,1)
            rcend1(ib1) = conformal_convert(non_commensurate,
     .       1,rcend(ib),rcpos(ib),nx,ny,nx1,ny1,3)
          else if (crcstra(ib).eq.'W') then
            rcpos1(ib1) = conformal_convert(non_commensurate,
     .       1,rcpos(ib),rcstart(ib),nx,ny,nx1,ny1,3)
            rcstart1(ib1) = conformal_convert(non_commensurate,
     .       2,rcpos(ib),rcstart(ib),nx,ny,nx1,ny1,0)
            rcend1(ib1) = conformal_convert(non_commensurate,
     .       2,rcpos(ib),rcend(ib),nx,ny,nx1,ny1,2)
          else if (crcstra(ib).eq.'E') then
            rcpos1(ib1) = conformal_convert(non_commensurate,
     .       1,rcpos(ib),rcstart(ib),nx,ny,nx1,ny1,1)
            rcstart1(ib1) = conformal_convert(non_commensurate,
     .       2,rcpos(ib),rcstart(ib),nx,ny,nx1,ny1,0)
            rcend1(ib1) = conformal_convert(non_commensurate,
     .       2,rcpos(ib),rcend(ib),nx,ny,nx1,ny1,2) 
          else if (crcstra(ib).eq.'V' .or. crcstra(ib).eq.'T' .or. 
     .             crcstra(ib).eq.'U' .or. crcstra(ib).eq.'C') then
            rcpos1(ib1) = 0
            rcstart1(ib1) = 0
            rcend1(ib1) = 0
          endif
          species_start1(ib1)=
     .     is1beginrange(species_start(ib),species_end(ib),ns,ns1,is1)
          species_end1(ib1)=
     .     is1endrange(species_start(ib),species_end(ib),ns,ns1,is1)
          cftmp(0:ns-1) = recyc(0:ns-1,ib)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          recyc1(0:ns1-1,ib1)=cftmp(0:ns1-1)
          cftmp(0:ns-1) = erecyc(0:ns-1,ib)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          erecyc1(0:ns1-1,ib1)=cftmp(0:ns1-1)
          cftmp(0:ns-1) = rcion(0:ns-1,ib)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          rcion1(0:ns1-1,ib1)=cftmp(0:ns1-1)
          cftmp(0:ns-1) = phys_sput(0:ns-1,ib)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          phys_sput1(0:ns1-1,ib1)=cftmp(0:ns1-1)
          cftmp(0:ns-1) = chem_sput(0:ns-1,ib)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          chem_sput1(0:ns1-1,ib1)=cftmp(0:ns1-1)
          chemsp1(ib1)=chemsp(ib)
          recyceir1(ib1)=recyceir(ib)
          volrecstart1(ib1)=volrecstart(ib)
          userfluxparm1(ib1,1:2)=userfluxparm(ib,1:2)
         endif
        enddo
        nstrai1 = ib1
        call convert_bundle_i
     .   (ns1,ns,nsm,is_neutral1,bundle_frac,track_index)
        do is = 0, ns1-1
          if (track_index(is).eq.itrack1) then
            track_index(is) = 1
          elseif (track_index(is).eq.itrack2) then
            track_index(is) = 2
          elseif (track_index(is).eq.itrack3) then
            track_index(is) = 3
          else
            track_index(is) = 0
          endif
        end do
        inquire(file='input.dat',exist=ex)
        if (.not.ex) then  ! must fill new arrays b2espcr1 and eb2spcr1 not provided by modified input.dat file
          iss = 0
          do is = 1, nspecies
            if (is1(eb2spcr(is)).ne.-1) then
              iss = iss + 1
              eb2spcr1(iss) = is1(eb2spcr(is))
              b2espcr1(
     .   is1beginrange(eb2spcr(is),eb2spcr(is)+nfluids(is),ns,ns1,is1):
     .     is1endrange(eb2spcr(is),eb2spcr(is)+nfluids(is),ns,ns1,is1))=
     .         is   
            end if
          end do
        else if (.not.same_bundling) then
          iss = 0
          do ib = 0, ns1-1
            is = 0
            found = .false.
            do while (.not.found .and. is.lt.ns)
              found = ib.eq.is1(is)
              is = is + 1
            end do
            if (.not.found) then
              if (is_neutral1(ib)) iss = iss+1
              b2espcr1(ib) = iss
              if (is_neutral1(ib)) eb2spcr1(iss) = ib
              if (.not.is_neutral1(ib) .and. nint(zn1(ib)).ne.6) then
                RECYC1(ib,:)=1.0_R8
                ERECYC1(ib,:)=0.3_R8
              endif
            endif
          end do
        endif
        do is = 1, nspecies
          if (eb2spcr1(is).ne.-1) then
            isp1(is) = b2espcr1(eb2spcr1(is))
          else
            isp1(is) = 0
          end if
        end do
        mlcmp1 = 0.0_R8
        do is = 1, nspecies
          if (isp1(is).ne.0) then
            mlcmp1(isp1(is),1:nmol) = mlcmp(is,1:nmol)
          endif
        end do
        targsp1 = 0
        do is = 1, ntrack
          ib1 = 0
          do ib = 1, nstrai
            if (savestratum(ib).or.crcstra(ib).eq.'T') then
              ib1 = ib1+1
              if (targsp(ib,is).ge.0) targsp1(ib1,is)=is1(targsp(ib,is))
            endif
          enddo
        enddo
        neutrals_filename = 
     .   '../'//neutrals_filename(1:lnblnk(neutrals_filename))//'.2'
        write(*,*) 'Writing ',trim(neutrals_filename)
        open(34,file=neutrals_filename)
        write(34,'(a)') '&NEUTRALS'
        write(34,'(a,i2,a)') 'NSTRAI=', nstrai1,','
        if (eirene_step_cpu.gt.0.0_R8)
     >   write(34,'(a,1pg9.2,a)') 'EIRENE_STEP_CPU=',eirene_step_cpu,','
        write(34,'(a,1pg9.2,a)') 'EIRENE_STEP_DT= ',eirene_step_dt,','
        if (eirene_mod.gt.1)
     >   write(34,'(a,i2,a)') 'EIRENE_MOD=', eirene_mod,','
        line = 'CRCSTRA=      '
        len=14
        do ib = 1, nstrai1
          line = line(1:len)//'  '''//crcstra1(ib)//''','
          len=len+3+1+2
        end do
        write(34,'(a)') line
        write(hlp_frm,'(a,i2,a)') '(a,',nstrai1,'(1x,i4,'',''))'
        write(34, hlp_frm) 'RCPOS=        ',(rcpos1(ib)  ,ib=1,nstrai1)
        write(34, hlp_frm) 'RCSTART=      ',(rcstart1(ib),ib=1,nstrai1)
        write(34, hlp_frm) 'RCEND=        ',(rcend1(ib)  ,ib=1,nstrai1)
        if (maxval(species_start1(1:nstrai1)).gt.0)
     >   write(34, hlp_frm) 
     .             'SPECIES_START=',(species_start1(ib),ib=1,nstrai1)
        if (minval(species_end1(1:nstrai1)).lt.ns1-1)
     >   write(34, hlp_frm) 
     .             'SPECIES_END=  ',(species_end1(ib),ib=1,nstrai1)
        write(hlp_frm,'(a,i2,a)') '(a,i2,a,1p,',ns1,'(1x,g9.2,'',''))'
        do ib = 1, nstrai1
          write(34, hlp_frm) 
     .     'RECYC(0,',ib,')= ',(recyc1(is,ib),is=0,ns1-1)
          write(34, hlp_frm) 
     .     'ERECYC(0,',ib,')=',(erecyc1(is,ib),is=0,ns1-1)
          if (maxval(rcion1(0:ns1-1,ib)).gt.0.0_R8)
     >     write(34, hlp_frm) 
     .     'RCION(0,',ib,')= ',(rcion1(is,ib),is=0,ns1-1)
        end do
        write(hlp_frm,'(a,i2,a)') '(a,i2,a,1p,',ns1,'(1x,g9.2,'',''))'
        todo = .false.
        do ib = 1, nstrai1
          todo = todo .or. chemsp1(ib)
        enddo
        write(hlp_frm,'(a,i2,a)') '(a,',nstrai1,'l1,a)'
        if (todo) write(34,hlp_frm) 'CHEMSP=',chemsp1(1:nstrai1),','
        write(hlp_frm,'(a,i2,a)') '(a,',ns1,'(1x,i2,'',''))'
        write(34, hlp_frm) 'TRACK_INDEX=',(track_index(is),is=0,ns1-1)
        track_chem_sput1 = .false.
        if (itrack1.gt.0) track_chem_sput1(1) = track_chem_sput(itrack1)
        if (itrack2.gt.0) track_chem_sput1(2) = track_chem_sput(itrack2)
        if (itrack3.gt.0) track_chem_sput1(3) = track_chem_sput(itrack3)
        write(hlp_frm,'(a,i2,a)') '(a,',itrack,'l1,a)'
        if (itrack.gt.0) write(34, hlp_frm) 
     .   'TRACK_CHEM_SPUT=',track_chem_sput1(1:itrack),','
        write(hlp_frm,'(a,i2,a)') '(a,i2,a,',nstrai1,'(1x,i2,'',''))'
        do is = 1, ntrack
          if (maxval(targsp1(1:nstrai1,is)).gt.0) write(34, hlp_frm) 
     .     'TARGSP(1,',is,')=',(targsp1(ib,is),ib=1,nstrai1)
        end do
        write(hlp_frm,'(a,i2,a)') '(a,1p,',nstrai1,'(1x,g9.2,'',''))'
        if (minval(recyceir1(1:nstrai1)).lt.1.0_R8)
     >   write(34, hlp_frm) 'RECYCEIR=',(recyceir1(1:nstrai1))
        write(hlp_frm,'(a,i2,a)') '(a,i2,a,1p,',ns1,'(1x,g9.2,'',''))'
        do ib = 1, nstrai1
          if (maxval(phys_sput1(0:ns1-1,ib)).gt.0.0_R8)
     >     write(34, hlp_frm) 
     .     'PHYS_SPUT(0,',ib,')=',(phys_sput1(is,ib),is=0,ns1-1)
          if (maxval(chem_sput1(0:ns1-1,ib)).gt.0.0_R8)
     >     write(34, hlp_frm) 
     .     'CHEM_SPUT(0,',ib,')=',(chem_sput1(is,ib),is=0,ns1-1)
        end do
        if (chemical_erosion_redep_fac.ne.1.0_R8)
     >   write(34, '(a,1pg9.2,a)')
     .    'CHEMICAL_EROSION_REDEP_FAC=',
     .     chemical_erosion_redep_fac,','
        if (chemical_erosion_be_fac) then
          write(34, '(a)')
     >     'CHEMICAL_EROSION_BE_FAC=T,'
          if (chemical_erosion_be_fac_a.ne.0.0_R8)
     >     write(34, '(a,1pg9.2,a)')
     .     'CHEMICAL_EROSION_BE_FAC_A=',
     .      chemical_erosion_be_fac_a,','
          if (chemical_erosion_be_fac_b.ne.0.0_R8)
     >     write(34, '(a,1pg9.2,a)')
     .     'CHEMICAL_EROSION_BE_FAC_B=',
     .      chemical_erosion_be_fac_b,','
          if (chemical_erosion_be_fac_c.ne.0.0_R8)
     >     write(34, '(a,1pg9.2,a)')
     .     'CHEMICAL_EROSION_BE_FAC_C=',
     .      chemical_erosion_be_fac_c,','
        endif
        write(hlp_frm,'(a,i2,a)') '(a,1p,',nstrai1,'(1x,g9.2,'',''))'
        do ib = 1, nstrai1
         if (crcstra1(ib).eq.'V' .and. volrecstart1(ib).ne.1.0e21_R8)
     >    write(34, '(a,i2,a,1pg9.2,a)') 
     .    'VOLRECSTART(',ib,')=',volrecstart1(ib),','
        end do
        if (volrecinc.ne.0.0_R8) write(34, '(a,1pg9.2,a)') 
     .   'VOLRECINC=',volrecinc,','
        write(34, '(a,1pg9.2,a)') 
     .   'VOLRECWT=',volrecwt,','
        do ib = 1, 4
          if (rf_neut(ib).ne.1.0_R8) write(34,'(a,i1,a,1pg9.2,a)') 
     .     'RF_NEUT(',ib,')=',rf_neut(ib)
        end do
        do ib = 1, nstrai1
          if (userfluxparm1(ib,1).ne.0.0_R8) 
     >     write(34, '(a,i2,a,1pg9.2,a)') 
     .     'USERFLUXPARM(',ib,',1)=',userfluxparm1(ib,1),','
          if (userfluxparm1(ib,2).ne.0.0_R8) 
     >     write(34, '(a,i2,a,1pg9.2,a)') 
     .     'USERFLUXPARM(',ib,',2)=',userfluxparm1(ib,2),','
        end do
        ib1 = 0
        do ib = 1, nstrai
          if (savestratum(ib).or.crcstra(ib).eq.'T') then
            ib1 = ib1+1
            do is = 1, nspecies
              if (isp1(is).ne.0 .and. gpfc(is,ib).ne.0.0_R8) 
     .         write(34,'(a,i1,a,i2,a,f5.2,a)')
     .          'GPFC(',isp1(is),',',ib,')=',gpfc(is,ib),','
            end do
          end if
        end do
#ifndef NO_B2_CHEM_SPUT
        do ib = 0, nlimps
          if (chemical_sputter_yield(ib).ne.0.0_R8)
     >     write(34,'(a,i3,a,1pg9.2,a)') 
     .     'CHEMICAL_SPUTTER_YIELD(',ib,')=',
     .     chemical_sputter_yield(ib),','
        end do
        if (fchar_chemical.ne.0.0_R8) write(34, '(a,1pg9.2,a)')
     .   'FCHAR_CHEMICAL=',fchar_chemical,','
        if (igass_chemical.ne.0) write(34, '(a,i2,a)')
     .   'IGASS_CHEMICAL=',igass_chemical,','
        if (itsput_chemical.ne.0) write(34, '(a,i2,a)')
     .   'ITSPUT_CHEMICAL=',itsput_chemical,','
        if (issput_chemical.ne.0) write(34, '(a,i2,a)')
     .   'ISSPUT_CHEMICAL=',issput_chemical,','
#endif
        if (neutrals_time_mod.ne.0.0_R8) write(34,'(a,f10.7,a)') 
     >   'NEUTRALS_TIME_MOD=', neutrals_time_mod,','
        if (neutrals_time_switch.ne.0.0_R8) write(34,'(a,f10.7,a)') 
     >   'NEUTRALS_TIME_SWITCH=', neutrals_time_switch,','
        write(hlp_frm,'(a,i2,a)') '(a,i1,a,',nstrai1,'(1x,i2,'',''))'
        if (maxval(lstrascl).gt.0) then
          do is = 1, nspecies
            if (isp1(is).ne.0) then
              ib1=0
              do ib = 1,nstrai
                if (savestratum(ib).or.crcstra(ib).eq.'T') then
                  ib1 = ib1+1
                  lstrascl1(ib1,isp1(is))=lstrascl(ib,is)
                endif
              enddo
            endif
          enddo         
          do is = 1, nspecies
            if (isp1(is).ne.0) then
              write(34, hlp_frm) 'LSTRASCL(1,',isp1(is),')=',
     .         (lstrascl1(ib,isp1(is)),ib=1,nstrai1)
            end if
          end do
        endif
        write(hlp_frm,'(a,i2,a)') '(a,',nmol,'(1x,i2,'',''))'
        if (maxval(lmolscl).gt.0) 
     >   write(34,hlp_frm) 'LMOLSCL=',(lmolscl(is),is=1,nmol)
        write(hlp_frm,'(a,i2,a)') '(a,',ns1,'(1x,i2,'',''))'
        write(34,hlp_frm) 'B2ESPCR=',(b2espcr1(is),is=0,ns1-1)
        write(hlp_frm,'(a,i2,a)') '(a,',nspecies1,'(1x,i2,'',''))'
        write(34,hlp_frm) 'EB2SPCR=',(eb2spcr1(is),is=1,nspecies1)
        write(hlp_frm,'(a,i2,a)') '(a,i2,a,',nspecies1,'(1x,i2,'',''))'
        do ib = 1, nmol
          if (maxval(mlcmp1(:,ib)).gt.0) then
            write(34,hlp_frm) 
     .       'MLCMP(1,',ib,')=',(mlcmp1(is,ib),is=1,nspecies1)
          endif
        end do
        ii=minloc(lcns)
        if(lcns(ii(1)).gt.0) then !{
          ncns=nstsid
        else !}{
          ncns=ii(1)-1
        end if !}
        ii=minloc(ltns)
        if(ltns(ii(1)).gt.0) then !{
          ntns=nstsid
        else !}{
          ntns=ii(1)-1
        end if !}
        if (ncns.gt.0) then
          write(hlp_frm,'(a,i2,a)') '(a,',ncns,'(1x,i2,'',''))'
          write(34,hlp_frm) 'LCNS=',(lcns(ib),ib=1,ncns)
        endif
        if (ntns.gt.0) then
          write(hlp_frm,'(a,i2,a)') '(a,',ntns,'(1x,i2,'',''))'
          write(34,hlp_frm) 'LTNS=',(ltns(ib),ib=1,ntns)
        endif
        nsns = 0
        do ib = 1, nstrai
          if (crcstra(ib).eq.'N'.or. crcstra(ib).eq.'W'.or.
     &        crcstra(ib).eq.'E') nsns = nsns + 1
          if (crcstra(ib).eq.'S') then
            if (mod(region1(rcstart1(ib),rcpos1(ib),0),4).ne.1) then
              nsns = nsns + 1
            end if
          end if
        end do
        if (nsns.gt.0) then
          write(hlp_frm,'(a,i2,a)') '(a,',nsns,'(1x,i3,'',''))'
          if(maxval(lsns(1:nsns)).gt.0.or.minval(lsns(1:nsns)).lt.0)
     >     write(34,hlp_frm) 'LSNS=',(lsns(ib),ib=1,nsns)
          if(maxval(ksns(1:nsns)).gt.0.or.minval(ksns(1:nsns)).lt.0)
     >     write(34,hlp_frm) 'KSNS=',(ksns(ib),ib=1,nsns) 
        endif     
        write(34,'(a,1pg9.2,a)') 'NEUT_SCL_LIM=',neut_scl_lim,','
        if (n_spcsrf.gt.0) then
          write(34,'(a,i2,a)') 'N_SPCSRF=',n_spcsrf,','
          ii=minloc(l_spcsrf)
          write(hlp_frm,'(a,i3,a)') '(a,',ii(1),'(i3,'',''))'
          write(34,hlp_frm) 'L_SPCSRF=',(l_spcsrf(ib),ib=1,ii(1))
          line = 'SPS_ID=  '
          len=9
          do ib = 1, n_spcsrf
            line = line(1:len)//''''//sps_id(ib)//''','
            len=len+1+8+2
          end do
          write(34,'(a)') line
          write(hlp_frm,'(a,i2,a)') '(a,',n_spcsrf,'(i10,'',''))'
          write(34, hlp_frm) 'I_SPCSRF=',(i_spcsrf(ib),ib=1,n_spcsrf)
          write(34, hlp_frm) 'J_SPCSRF=',(j_spcsrf(ib),ib=1,n_spcsrf)
          write(hlp_frm,'(a,i2,a)') '(a,1p,',n_spcsrf,'(g10.3,'',''))'
          write(34, hlp_frm) 'SPS_ABSR=',(sps_absr(ib),ib=1,n_spcsrf)
          write(34, hlp_frm) 'SPS_TRNO=',(sps_trno(ib),ib=1,n_spcsrf)
          write(34, hlp_frm) 'SPS_TRNI=',(sps_trni(ib),ib=1,n_spcsrf)
          write(34, hlp_frm) 'SPS_MTRI=',(sps_mtri(ib),ib=1,n_spcsrf)
          write(34, hlp_frm) 'SPS_TMPR=',(sps_tmpr(ib),ib=1,n_spcsrf)
          write(34, hlp_frm) 'SPS_SPPH=',(sps_spph(ib),ib=1,n_spcsrf)
          write(34, hlp_frm) 'SPS_SPCH=',(sps_spch(ib),ib=1,n_spcsrf)
          write(hlp_frm,'(a,i2,a)') '(a,',n_spcsrf,'(i10,'',''))'
          write(34, hlp_frm) 'SPS_SGRP=',(sps_sgrp(ib),ib=1,n_spcsrf)
        endif
        do ib = 1, 100
          if (debug_flags(ib).ne.0) write(34,'(a,i3,a,i2,a)')
     .     'DEBUG_FLAGS(',ib,')=',debug_flags(ib),','
        end do
        if (dbg_eir_mc.ne.0) write(34,'(a,i2,a)') 
     .   'DBG_EIR_MC=',dbg_eir_mc,','
        write(34,'(a)') '/'
        close(34)
      endif      

*   ..we now do b2ah.dat, last because we need to overwrite zn, etc...
      call b2ahx0 (ninp(5), idum(0))
      call xertst (idum(0).eq.ns, 
     . 'incorrect species number in b2ah.dat')
      write(nout(3),'(a)') '*dimens         (ns; free format)'
      write(nout(3),*) ns1
      call b2rflb (ninp(5), lblcp(61:120))
      write(nout(3),'(a)') 
     . '*label          (lblcp, character*60; free format)' 
      write(nout(3),'(a,a,a)') '''',trim(lblcp(61:120)),''''
      call b2rfcp (ninp(5), ns)
      zamin(0:ns1-1) = zamin1(0:ns1-1) ! need to do this transfer
      zamax(0:ns1-1) = zamax1(0:ns1-1) ! so that the write in
      zn(0:ns1-1) = zn1(0:ns1-1)       ! b2wfcp is correct
      am(0:ns1-1) = am1(0:ns1-1)
      do ireg = 0, cbnrso+cbnrno+cbnrwe+cbnrea-1
        do i = 0, 7
          cftmp(0:ns-1) = cbsna(i,0:ns-1,ireg)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cbsna(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
          cftmp(0:ns-1) = cbsmo(i,0:ns-1,ireg)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cbsmo(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
          cftmp(0:ns-1) = cbshi(i,0:ns-1,ireg)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cbshi(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
          cftmp(0:ns-1) = cbrec(i,0:ns-1,ireg)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cbrec(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
          cftmp(0:ns-1) = cbmsb(i,0:ns-1,ireg)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cbmsb(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
        enddo
      enddo
      do i = 0, 7
        cftmp(0:ns-1) = cfdf0(i,0:ns-1)
        call convert_bundle_r
     .   (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
        cfdf0(i,0:ns1-1) = cftmp(0:ns1-1)
        cftmp(0:ns-1) = cfdna(i,0:ns-1)
        call convert_bundle_r
     .   (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
        cfdna(i,0:ns1-1) = cftmp(0:ns1-1)
        cftmp(0:ns-1) = cfdpa(i,0:ns-1)
        call convert_bundle_r
     .   (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
        cfdpa(i,0:ns1-1) = cftmp(0:ns1-1)
        cftmp(0:ns-1) = cfvla(i,0:ns-1)
        call convert_bundle_r
     .   (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
        cfvla(i,0:ns1-1) = cftmp(0:ns1-1)
        cftmp(0:ns-1) = cfvsa(i,0:ns-1)
        call convert_bundle_r
     .   (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
        cfvsa(i,0:ns1-1) = cftmp(0:ns1-1)
        cftmp(0:ns-1) = cfhci(i,0:ns-1)
        call convert_bundle_r
     .   (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
        cfhci(i,0:ns1-1) = cftmp(0:ns1-1)
      enddo
      call b2wfcp (nout(3), 0, ns1)

*   ..we now write a new SOLPS_FLUIDS file if needed
      if (.not.same_bundling) then
        write(*,*) 'Writing SOLPS_FLUIDS.2'
        allocate(species_name(0:ns1-1))
        do is = 0, ns1-1
          call species(is,species_name(is))
        enddo
        write(hlp_frm,'(a1,i2,a7)') '(',ns1,'(a,1x))'
        open(unit=99,file='../SOLPS_FLUIDS.2')
        write(99,hlp_frm) (trim(species_name(is)),is=0,ns1-1)
        close(99)
        deallocate(species_name)
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab
     1 ('error trying to read non-default internal parameters')
      return

*-----------------------------------------------------------------------
*.formats

  101 format (/'  The new grid contains nx1*ny1 interior cells, ',
     &  'where nx1 =',i4,' and ny1 =',i4,'.')
 1028 format ()
 6662 FORMAT (L1,1X,A24,1X,I1,1X,4(2I3,1X))
 6664 FORMAT (1P,6E12.4)
 6665 FORMAT (12(5L1,1X))
 6666 FORMAT (12I6)
66661 FORMAT (1P,I3,1X,A6,1X,A4,A9,A3,2I3,3E12.4)
66662 FORMAT (1P,2E12.4,10I6)
66664 FORMAT (1P,I6,6X,5E12.4)
66665 FORMAT (1P,2I6,4E12.4)
66666 FORMAT (I2,1X,A8,1X,12(I2,1X))

*-----------------------------------------------------------------------
*.end b2ytdr

      end subroutine b2ytdr

      integer function convert_position(i,nold,nnew,iconf,atend)
      implicit none
      integer i,nold,nnew,iconf
      logical atend
      intrinsic max

cxpb  This function converts a cell position i when the mesh 
cxpb  is changing from size nold to size nnew. 
cxpb  If atend is .true., then we want the right (or top) 
cxpb  position of the range of cells accessible.
cxpb  iconf : configuration index : 0 : boundary cells on both sides
cxpb                                1 : boundary cell on beginning side only
cxpb                                2 : boundary cell on end side only
cxpb                                3 : no boundary cells

      if (i.eq.-1 .and. iconf.le.1) then
        convert_position = -1
      else if (i.eq.0 .and. iconf.ge.2) then
        convert_position = 0
      else if (i.eq.nold .and. (iconf.eq.0 .or. iconf.eq.2)) then
        convert_position = nnew
      else if (i.eq.nold-1 .and. (iconf.eq.1 .or. iconf.eq.3)) then
        convert_position = nnew-1
      else if (.not.atend) then
        convert_position = i*nnew/nold
      else if (atend) then
        convert_position = max(0,(i+1)*nnew/nold-1)
      endif

      return
      end function convert_position

      integer function convert_eirene        
     . (non_commensurate,xory,ix,iy,nx,ny,nx1,ny1,iside)
      use b2mod_geo2
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer non_commensurate,xory,ix,iy,nx,ny,nx1,ny1,iside
c xory  = 1 : x-direction
c       = 2 : y-direction
c iside = 0 : bottom
c       = 1 : left
c       = 2 : top
c       = 3 : right
c       = 4 : center
*   ..local variables
      integer ixpos, iypos, icut
      integer eirleft(6), eirright(6)
      integer nxtl, nxtr
*   ..procedures
      integer conformal_convert
      external conformal_convert, get_nxt
cxpb Do not forget that in Eirene X is Radial and Y is poloidal !

      iypos = ix-1 ! Eirene IX is B2 IY
      eirleft = nx+7
      eirright= nx+7
      if (nncut.eq.2) call get_nxt(nx, nxtl, nxtr)
      if (non_commensurate.eq.1) then
        if (nncut.eq.0) then
          eirleft(1) = 1
          eirright(1)= nx+1
        else if (nncut.eq.1) then
          eirleft(1) = 1
          eirright(1)= leftcut(1)+1
          eirleft(2) = leftcut(1)+2
          eirright(2)= rightcut(1)+2
          eirleft(3) = rightcut(1)+3
          eirright(3)= nx+3
        else if (nncut.eq.2) then
          eirleft(1) = 1
          eirright(1)= leftcut(1)+1
          eirleft(2) = leftcut(1)+2
          eirright(2)= leftcut(2)+2
          eirleft(3) = leftcut(2)+3
          eirright(3)= nxtl+3
          eirleft(4) = nxtr+3
          eirright(4)= rightcut(2)+3
          eirleft(5) = rightcut(2)+4
          eirright(5)= rightcut(1)+4
          eirleft(6) = rightcut(1)+5
          eirright(6)= nx+5
        endif
      else if (non_commensurate.eq.0) then  
c-xpb  Only the new leftcut and rightcut geometry arrays are defined here
        if (nncut.eq.0) then
          eirleft(1) = 1
          eirright(1)= nx+1
        else if (nncut.eq.1) then
          eirleft(1) = 1
          eirright(1)= (leftcut1(1)*nx)/nx1+1
          eirleft(2) = (leftcut1(1)*nx)/nx1+2
          eirright(2)= (rightcut1(1)*nx)/nx1+2
          eirleft(3) = (rightcut1(1)*nx)/nx1+3
          eirright(3)= nx+3
        else if (nncut.eq.2) then
          eirleft(1) = 1
          eirright(1)= (leftcut1(1)*(nx-2))/(nx1-2)+1
          eirleft(2) = (leftcut1(1)*(nx-2))/(nx1-2)+2
          eirright(2)= (leftcut1(2)*(nx-2))/(nx1-2)+2
          eirleft(3) = (leftcut1(2)*(nx-2))/(nx1-2)+3
          eirright(3)= nxtl+3
          eirleft(4) = nxtr+3
          eirright(4)= (rightcut1(2)*(nx-2))/(nx1-2)+3
          eirleft(5) = (rightcut1(2)*(nx-2))/(nx1-2)+4
          eirright(5)= (rightcut1(1)*(nx-2))/(nx1-2)+4
          eirleft(6) = (rightcut1(1)*(nx-2))/(nx1-2)+5
          eirright(6)= nx+5
        endif
      else
        write(*,*) 
     .   'In convert_eirene = 0 : non_commensurate test failed'
        convert_eirene = 0
        return
      endif
      icut = 0
      do while (iy.ge.eirright(icut+1).and.icut.le.5)
        icut = icut+1
      enddo
      if (nncut.eq.2 .and. icut.gt.3) icut = icut - 1
      ixpos = iy-1-icut ! Eirene IY is B2 IX
      if (xory.eq.1) then
        convert_eirene = 1 + conformal_convert
     .   (non_commensurate,2,ixpos,iypos,nx,ny,nx1,ny1,iside)
      else if (xory.eq.2) then
        convert_eirene = icut + 1 + conformal_convert
     .   (non_commensurate,1,ixpos,iypos,nx,ny,nx1,ny1,iside)
      else
        write(*,*) 'In convert_eirene = 0 : xory test failed'
        convert_eirene = 0
      endif

      return
      end function convert_eirene

      integer function conformal_convert
     . (non_commensurate,xory,ix,iy,nx,ny,nx1,ny1,iside)
      use b2mod_indirect
      use b2mod_geo2
      implicit none
*   ..input arguments (unchanged on exit)
      integer non_commensurate,xory,ix,iy,nx,ny,nx1,ny1,iside
*   ..local variables
      integer ireg, ipos, nsize, nsize1
      integer ixiso, ixiso1, inseliy, inseliy1
      logical atend
*   ..procedures
      integer convert_position
      external convert_position
      intrinsic min, max

c xory  = 1 : x-direction
c       = 2 : y-direction
c iside = 0 : bottom
c       = 1 : left
c       = 2 : top
c       = 3 : right
c       = 4 : center

      atend = (xory.eq.1 .and. iside.eq.3) .or. 
     &        (xory.eq.2 .and. iside.eq.2)
      if (xory.eq.1) then
       ipos = ix
       nsize = nx
       nsize1 = nx1
      else if (xory.eq.2) then
       ipos = iy
       nsize = ny
       nsize1 = ny1
      else
       conformal_convert=-2
       return
      endif
      if (non_commensurate.eq.0) then
       conformal_convert=convert_position(ipos,nsize,nsize1,0,atend)
      else
       if (ipos.eq.-1) then  ! do easy cases
        conformal_convert=-1
       else if (ipos.eq.nsize) then
        conformal_convert=nsize1
       else      
        ireg=region(ix,iy,0) 
        if (nncut.eq.2) then  ! look for isolating cut
         if (ireg.eq.4) then
          ixiso = leftcut(2)
          do while (region(ixiso+1,topcut(2),0).eq.ireg)
           ixiso = ixiso+1
          end do
          ixiso1 = leftcut1(2)
          do while (region1(ixiso1+1,topcut1(2),0).eq.ireg)
           ixiso1 = ixiso1+1
          end do
         else if (ireg.eq.7) then
          ixiso = rightcut(2)-1
          do while (region(ixiso-1,topcut(2),0).eq.ireg)
           ixiso = ixiso-1
          end do
          ixiso1 = rightcut1(2)
          do while (region1(ixiso1-1,topcut1(2),0).eq.ireg)
           ixiso1 = ixiso1-1
          end do
         else
          ixiso = -2
          ixiso1 = -2
         endif
        else
         ixiso = -2
         ixiso1 = -2
        endif 
        if (periodic_bc.eq.1 .and. nncut.eq.1) then    ! look for inseliy
         if (ireg.eq.2) then
          inseliy = topcut(1)+1
          do while (region(nx/2,inseliy+1,0).eq.ireg)
           inseliy = inseliy+1
          enddo
          inseliy1 = topcut1(1)+1
          do while (region1(nx1/2,inseliy1+1,0).eq.ireg)
           inseliy1 = inseliy1+1
          enddo
         else if (ireg.eq.5) then
          inseliy = ny
          do while (region(nx/2,inseliy-1,0).eq.ireg)
           inseliy = inseliy-1
          enddo
          inseliy1 = ny1
          do while (region1(nx1/2,inseliy1-1,0).eq.ireg)
           inseliy1 = inseliy1-1
          enddo
         else
          inseliy = -2
          inseliy1 = -2
         endif
        else
         inseliy = -2
         inseliy1 = -2
        endif
        if (xory.eq.1) then   ! do other easy ix cases
         if (ix.eq.0) then
          conformal_convert=0
         else if (ix.eq.leftcut(1)-1) then
          conformal_convert=leftcut1(1)-1
         else if (ix.eq.leftcut(1)) then
          conformal_convert=leftcut1(1)
         else if (ix.eq.leftcut(2)-1) then
          conformal_convert=leftcut1(2)
         else if (ix.eq.leftcut(2)-1) then
          conformal_convert=leftcut1(2)
         else if (ix.eq.rightcut(1)-1) then
          conformal_convert=rightcut1(1)-1
         else if (ix.eq.rightcut(1)) then
          conformal_convert=rightcut1(1)
         else if (ix.eq.rightcut(2)-1) then
          conformal_convert=rightcut1(2)
         else if (ix.eq.rightcut(2)-1) then
          conformal_convert=rightcut1(2)
         else if (ix.eq.nx-1) then
          conformal_convert=nx1-1
         else if (ix.eq.ixiso-1) then
          conformal_convert=ixiso1-1
         else if (ix.eq.ixiso) then
          conformal_convert=ixiso1
         else if (ix.eq.ixiso+1) then
          conformal_convert=ixiso1+1
         else if (nncut.eq.0 .and. periodic_bc.eq.0) then  ! 2-D slab case
          nsize = nx
          nsize1 = nx1
          conformal_convert = convert_position(ix,nsize,nsize1,0,atend)
         else if (nncut.eq.0 .and. ireg.eq.1) then    ! limiter core region
          nsize = nx
          nsize1 = nx1
          conformal_convert = convert_position(ix,nsize,nsize1,3,atend)
         else if (periodic_bc.eq.1 .and. ireg.eq.2) then  ! limiter or island SOL
          nsize = nx
          nsize1 = nx1
          conformal_convert = convert_position(ix,nsize,nsize1,0,atend)
         else if (periodic_bc.eq.1 .and. ireg.eq.5) then  ! island
          nsize = nx
          nsize1 = nx1
          conformal_convert = convert_position(ix,nsize,nsize1,3,atend)
         else if (ix.gt.0 .and. ix.lt.leftcut(1)-1) then
          nsize = leftcut(1)
          nsize1 = leftcut1(1)
          conformal_convert = convert_position(ix,nsize,nsize1,1,atend)
         else if (nncut.eq.1 .and. 
     .            ix.gt.leftcut(1) .and. ix.lt.rightcut(1)-1) then
          nsize = rightcut(1)-leftcut(1)
          nsize1 = rightcut1(1)-leftcut1(1)
          conformal_convert = leftcut1(1) + 
     .     convert_position(ix-leftcut(1),nsize,nsize1,3,atend)
         else if (nncut.eq.2 .and. 
     .            ix.gt.leftcut(1) .and. ix.le.leftcut(2)-1) then
          nsize = leftcut(2)-leftcut(1)
          nsize1 = leftcut1(1)-leftcut1(1)
          conformal_convert = leftcut1(1) + 
     .     convert_position(ix-leftcut(1),nsize,nsize1,3,atend)
         else if (nncut.eq.2 .and. 
     .            ix.gt.leftcut(2) .and. ix.lt.ixiso) then
          nsize = ixiso-leftcut(2)
          nsize1 = ixiso1-leftcut1(2)
          conformal_convert = leftcut1(2) +
     .     convert_position(ix-leftcut(2),nsize,nsize1,2,atend)
         else if (nncut.eq.2 .and.
     .            ix.gt.ixiso .and. ix.lt.rightcut(2)-1) then
          nsize = rightcut(2)-ixiso-1
          nsize1 = rightcut1(2)-ixiso1-1
          conformal_convert = ixiso1 +
     .     convert_position(ix-ixiso-1,nsize,nsize1,1,atend)
         else if (nncut.eq.2 .and. 
     .            ix.gt.rightcut(2) .and. ix.le.rightcut(1)-1) then
          nsize = rightcut(1)-rightcut(2)
          nsize1 = rightcut1(1)-rightcut1(2)
          conformal_convert = rightcut1(2) +
     .     convert_position(ix-rightcut(2),nsize,nsize1,3,atend)
         else if (ix.gt.rightcut(1) .and. ix.lt.nx-1) then
          nsize = nx-rightcut(1)-1
          nsize1 = nx1-rightcut1(1)-1
          conformal_convert = rightcut1(1) +
     .     convert_position(ix-rightcut(1),nsize,nsize1,2,atend)
         else
          write(*,*) 'ix conformal conversion not succesful !'
          conformal_convert = -2
         endif
        else if (xory.eq.2) then   ! do other easy iy cases
         if (iy.eq.0) then
          conformal_convert=0
         else if (nncut.ge.1 .and. iy.eq.topcut(1)) then
          conformal_convert=topcut1(1)
         else if (nncut.ge.1 .and. iy.eq.topcut(1)+1) then
          conformal_convert=topcut1(1)+1
         else if (nncut.eq.2 .and. iy.eq.topcut(2)) then
          conformal_convert=topcut1(2)
         else if (nncut.eq.2 .and. iy.eq.topcut(2)+1) then
          conformal_convert=topcut1(2)+1
         else if (nncut.eq.1 .and. periodic_bc.eq.1 .and.
     .    iy.eq.inseliy) then
          conformal_convert=inseliy1
         else if (iy.eq.ny-1) then
          conformal_convert=ny1-1
         else if (nncut.eq.0) then
          nsize = ny
          nsize1 = ny1
          conformal_convert = convert_position(iy,nsize,nsize1,0,atend)
         else if (nncut.eq.1 .and. iy.lt.topcut(1)) then
          nsize = topcut(1)
          nsize1 = topcut1(1)
          conformal_convert = convert_position(iy,nsize,nsize1,2,atend)
         else if (nncut.eq.2 .and. iy.lt.min(topcut(1),topcut(2))) then
          nsize = min(topcut(1),topcut(2))
          nsize1 = min(topcut1(1),topcut1(2))
          conformal_convert = convert_position(iy,nsize,nsize1,2,atend)
         else if (nncut.eq.2 .and. iy.gt.min(topcut(1),topcut(2)) .and.
     .                             iy.lt.max(topcut(1),topcut(2))) then
          nsize = max(topcut(1),topcut(2))-min(topcut(1),topcut(2))
          nsize1 = max(topcut1(1),topcut1(2))-min(topcut1(1),topcut1(2))
          conformal_convert = min(topcut1(1),topcut1(2)) + 
     .     convert_position(iy-min(topcut(1),topcut(2)),
     .      nsize,nsize1,3,atend)
         else if (nncut.eq.2 .and. iy.gt.max(topcut(1),topcut(2))) then
          nsize = ny-max(topcut(1),topcut(2))
          nsize1 = ny1-max(topcut1(1),topcut1(2))
          conformal_convert = max(topcut1(1),topcut1(2)) + 
     .     convert_position(iy-max(topcut(1),topcut(2)),
     .      nsize,nsize1,2,atend)
         else if (nncut.eq.1 .and. periodic_bc.eq.0 .and. 
     .            iy.gt.topcut(1)) then
          nsize = ny-topcut(1)
          nsize1 = ny1-topcut1(1)
          conformal_convert = topcut1(1) +
     .     convert_position(iy-topcut(1),nsize,nsize1,2,atend)
         else if (nncut.eq.1 .and. periodic_bc.eq.1 .and. 
     .            iy.gt.topcut(1) .and. iy.lt.inseliy) then
          nsize = inseliy-topcut(1)+1
          nsize1 = inseliy1-topcut1(1)+1
          conformal_convert = topcut1(1) +
     .     convert_position(iy-topcut(1),nsize,nsize1,3,atend)
         else if (nncut.eq.1 .and. periodic_bc.eq.1 .and. 
     .            iy.gt.inseliy) then
          nsize = ny-inseliy
          nsize1 = ny1-inseliy1
          conformal_convert = inseliy1 +
     .     convert_position(iy-inseliy,nsize,nsize1,2,atend)
         else
          write(*,*) 'iy conformal conversion not succesful !'
          conformal_convert = -2 
         endif
        endif
       endif 
      endif
      
      return
      end function conformal_convert

      subroutine de_bundle(ns, zn, am, zamin, zamax,
     *   is0, bundle_contr, zn0, am0, zamin0, zamax0)
      use b2mod_types
      implicit none
*   ..arguments for input bundles (unchanged upon exit)
      integer ns
      real (kind=R8) :: 
     *   zn(0:ns-1), am(0:ns-1), zamin(0:ns-1), zamax(0:ns-1)
*   ..arguments for de-bundled data (unspecified upon entry)
      integer is0(0:100)
      real (kind=R8) :: 
     *   zn0(0:100), am0(0:100), zamin0(0:100), zamax0(0:100),
     *   bundle_contr(0:ns-1,0:100)
*   ..local variables
      integer i, is, ic
*   ..procedures
      intrinsic nint
      external xertst

      ic = -1
      is0 = -1
      do is = 0, ns-1
        do i = nint(zamin(is)), nint(zamax(is))
          ic = ic + 1
          call xertst (ic.le.100, 'increase size of is0')
          is0(ic) = is
          zn0(ic) = zn(is)
          am0(ic) = am(is)
          zamin0(ic) = i
          zamax0(ic) = i
          bundle_contr(is,ic) = 1.0_R8 / (zamax(is)-zamin(is)+1.0_R8)
        end do 
      end do
cxpb Old bundle 'is' receives contribution 'bundle_contr(is,ic)' from the
cxpb basic species 'ic'

      return
      end subroutine de_bundle

      subroutine re_bundle (
     *   ns, is0, bundle_contr, zn0, am0, zamin0, zamax0,
     *   ns1, is1, bundle_frac, zn1, am1, zamin1, zamax1)
      use b2mod_types
      implicit none
*   ..arguments for input de-bundled data (unchanged upon exit)
      integer ns
      integer is0(0:100)
      real (kind=R8) :: 
     *   zn0(0:100), am0(0:100), zamin0(0:100), zamax0(0:100),
     *   bundle_contr(0:ns-1,0:100)
*   ..arguments for re-bundled data (unchanged upon entry)
      integer ns1
      real (kind=R8) :: 
     *   zn1(0:ns1-1), zamin1(0:ns1-1), zamax1(0:ns1-1)
*   ..re-bundled data to compute (unspecified upon entry)
      integer is1(0:ns-1)
      real (kind=R8) :: 
     *   am1(0:ns1-1), bundle_frac(0:ns1-1,0:ns-1)
*   ..local variables
      integer ic, iss, nz, ii(1)
      logical found
*   ..procedures
      intrinsic nint
      external xerrab

cxpb We first need to find the total number of basic species nz
      ii = minloc(is0)
      nz = ii(1)-1
      do ic = 0, nz-1
        found = .false.
        iss = 0
        do while (.not.found .and. iss.lt.ns1)
          found = nint(zn1(iss)).eq.nint(zn0(ic)) .and.
     &     nint(zamin1(iss)).le.nint(zamin0(ic)) .and.
     &     nint(zamax1(iss)).ge.nint(zamax0(ic))
          if (.not.found) iss = iss + 1
        enddo
        if (iss.lt.ns1) then
          is1(is0(ic)) = iss
          am1(iss) = am0(ic)
          bundle_frac(iss,is0(ic)) = bundle_frac(iss,is0(ic)) + 
     &     bundle_contr(is0(ic),ic)
cxpb New bundle 'iss' receives contribution 'bundle_frac(iss,is)'
cxpb from old bundle 'is'
        else
          is1(is0(ic)) = -1
          write(*,*) 'Species ',is0(ic),' will be removed !'
cxpb          call xerrab ('Bundle species not found !')
        endif
      enddo

      return
      end subroutine re_bundle

      subroutine convert_bundle_r
     . (ns1,ns,nsm,is_neutral1,bundle_frac,work)
      use b2mod_types
      use b2mod_b2cmpa
      implicit none
      integer ns1, ns, nsm
      integer is1, is, is0
      logical is_neutral1(0:ns1-1), same
      real (kind=R8) ::
     *   bundle_frac(0:ns1-1,0:ns-1), work(0:nsm-1),
     *   rold(0:ns-1), sum, value
      logical filled

cxpb We assign to the bundle the value from the highest species from that bundle

      rold(0:ns-1)=work(0:ns-1)
      work(0:ns1-1)=0.0_R8
      do is1 = 0, ns1-1
        filled = .false.
        sum = 0.0_R8
        do is = ns-1, 0, -1
          if (bundle_frac(is1,is).gt.0.0_R8) then
            sum = sum + bundle_frac(is1,is)
            if (.not.filled) then
              work(is1) = rold(is)
              filled = .true.
            endif
          endif
        enddo
        if (sum.eq.0.0_R8) then  ! this is a new species 
*  When all species have the same setting, we copy it unto new species
          if (maxval(rold).eq.minval(rold)) then
            work(is1)=maxval(rold)
          else
*  We check if all neutral or all ionised values are the same
*  and if so use that
            is0 = 0
            do while (is0.lt.ns .and. 
     &                is_neutral(is0).neqv.is_neutral1(is1))
              is0 = is0 + 1
            end do
            if (is0.lt.ns) then
              value = rold(is0)
              same = .true.
            else
              same = .false.
            endif
            do is = is0+1, ns-1
              if (is_neutral(is).eqv.is_neutral1(is1)) then
                same = same.and.value.eq.rold(is)
              endif
            enddo
            if (same) work(is1) = value
          endif
        endif
      enddo

      return
      end subroutine convert_bundle_r

      subroutine convert_bundle_i
     . (ns1,ns,nsm,is_neutral1,bundle_frac,iwrk)
      use b2mod_types
      use b2mod_b2cmpa
      implicit none
      integer ns1, ns, nsm, iwrk(0:nsm-1), iold(0:ns-1)
      integer is1, is, is0, ival
      logical is_neutral1(0:ns1-1)
      real (kind=R8) ::
     *   bundle_frac(0:ns1-1,0:ns-1), sum
      logical filled, same

cxpb We assign to the bundle the value from the highest species from that bundle

      iold(0:ns-1)=iwrk(0:ns-1)
      iwrk(0:ns1-1)=0
      do is1 = 0, ns1-1
        filled = .false.
        sum = 0.0_R8
        do is = ns-1, 0, -1
          if (bundle_frac(is1,is).gt.0.0_R8) then
            sum = sum + bundle_frac(is1,is)
            if (.not.filled) then
              iwrk(is1) = iold(is)
              filled = .true.
            endif
          endif
        enddo
        if (sum.eq.0.0_R8) then  ! this is a new species 
*  When all species have the same setting, we copy it unto new species
          if (maxval(iold).eq.minval(iold)) then
            iwrk(is1)=maxval(iold)
          else
*  We check if all neutral or all ionised values are the same
*  and if so use that
            is0 = 0
            do while (is0.lt.ns .and. 
     &                is_neutral(is0).neqv.is_neutral1(is1))
              is0 = is0 + 1
            end do
            if (is0.lt.ns) then
              ival = iold(is0)
              same = .true.
            else
              same = .false.
            endif
            do is = is0+1, ns-1
              if (is_neutral(is).eqv.is_neutral1(is1)) then
                same = same.and.ival.eq.iold(is)
              endif
            enddo
            if (same) iwrk(is1) = ival
          endif
        endif
      enddo

      return
      end subroutine convert_bundle_i

      subroutine convert_bundle_l
     . (ns1,ns,nsm,is_neutral1,bundle_frac,lwrk)
      use b2mod_types
      use b2mod_b2cmpa
      implicit none
      integer ns1, ns, nsm
      integer is1, is, is0
      logical is_neutral1(0:ns1-1)
      real (kind=R8) ::
     *   bundle_frac(0:ns1-1,0:ns-1), sum
      logical lwrk(0:nsm-1), lold(0:ns-1)
      logical filled, same, lval

cxpb We assign to the bundle the value from the highest species from that bundle

      lold(0:ns-1)=lwrk(0:ns-1)
      lwrk(0:ns1-1)=.true.
      do is1 = 0, ns1-1
        filled = .false.
        sum = 0.0_R8
        do is = ns-1, 0, -1
          if (bundle_frac(is1,is).gt.0.0_R8) then
            sum = sum + bundle_frac(is1,is)
            if (.not.filled) then
              lwrk(is1) = lold(is)
              filled = .true.
            endif
          endif
        enddo
        if (sum.eq.0.0_R8) then  ! this is a new species 
*  When all species have the same setting, we copy it unto new species
          same = .true.
          do is = 1, ns-1
            same = same .and. lold(0).eqv.lold(is)
          enddo
          if (same) lwrk(is1) = lold(0)
        else
*  We check if all neutral or all ionised values are the same
*  and if so use that
          is0 = 0
          do while (is0.lt.ns .and. 
     &              is_neutral(is0).neqv.is_neutral1(is1))
            is0 = is0 + 1
          end do
          if (is0.lt.ns) then
            lval = lold(is0)
            same = .true.
          else
            same = .false.
          endif
          do is = is0+1, ns-1
            if (is_neutral(is).eqv.is_neutral1(is1)) then
              same = same.and.lval.eqv.lold(is)
            endif
          enddo
          if (same) lwrk(is1) = lval
        endif
      enddo

      return
      end subroutine convert_bundle_l

      subroutine sum_bundle(ns1,ns,nsm,is_neutral1,bundle_frac,work)
      use b2mod_types
      use b2mod_b2cmpa
      implicit none
      integer ns1, ns, nsm
      integer is1, is, is0
      logical is_neutral1(0:ns1-1), same
      real (kind=R8) ::
     *   bundle_frac(0:ns1-1,0:ns-1), work(0:nsm-1),
     *   rold(0:ns-1), sum, value

cxpb We assign to the new bundle a weighted sum of the values from the old bundle

      rold(0:ns-1)=work(0:ns-1)
      work(0:ns1-1)=0.0_R8
      do is1 = 0, ns1-1
        sum = 0.0_R8
        do is = 0, ns-1
          sum = sum + bundle_frac(is1,is)
          work(is1) = work(is1) + bundle_frac(is1,is)*rold(is)
        enddo
        if (sum.eq.0.0_R8) then  ! this is a new species 
*  When all species have the same setting, we copy it unto new species
          if (maxval(rold).eq.minval(rold)) then
            work(is1)=maxval(rold)
          else
*  We check if all neutral or all ionised values are the same
*  and if so use that
            is0 = 0
            do while (is0.lt.ns .and. 
     &                is_neutral(is0).neqv.is_neutral1(is1))
              is0 = is0 + 1
            end do
            if (is0.lt.ns) then
              value = rold(is0)
              same = .true.
            else
              same = .false.
            endif
            do is = is0+1, ns-1
              if (is_neutral(is).eqv.is_neutral1(is1)) then
                same = same.and.value.eq.rold(is)
              endif
            enddo
            if (same) work(is1) = value
          endif
        endif
      enddo

      return
      end subroutine sum_bundle

      integer function B2IDEZ(I,J,IZIF)
      implicit none
      integer i, j, izif
C
C   I IS AN INTEGER WITH N DECIMALS, THE FUNCTION RETURNS
C   IDEZ, THE J-TH DECIMAL
C   J MUST BE .LE. IZIF
C      IF N .LT. IZIF AND J .GT. N  , IDEZ=0
C      IF N .GT. IZIF AND J .EQ. IZIF  , IDEZ=PART OF I ON THE LEFT OF
C      THE DECIMAL J, INCLUDING THE J-TH DECIMAL
C      E.G.: I=1234,IZIF=3, J=1: IDEZ=4
C                           J=2: IDEZ=3
C                           J=3: IDEZ=12
C                           J=4: ERROR, J GT IZIF
C
      integer idif, iz, k, it, iq
      external subini, subend, xertst

      call subini('b2idez')
      IDIF=IZIF-J+1
      call xertst(idif.gt.0, 'ERROR IN FUNCTION B2IDEZ') 
      IZ=I
      DO 1 K=1,IDIF
        IT=10**(IZIF-K)
        IQ=IZ/IT
        IZ=IZ-IQ*IT
1     CONTINUE
      B2IDEZ=IQ
      call subend()
      return
      end function b2idez

      integer function is1beginrange(isbegin,isend,ns,ns1,is1)
      implicit none
      integer isbegin, isend, ns, ns1
      integer is1(0:ns-1)
      integer is
      intrinsic min

      if (isbegin.eq.0 .and. isend.eq.ns-1) then
        is1beginrange = 0
      else if (isbegin.eq.isend .and. is1(isbegin).ne.-1) then
        is1beginrange = is1(isbegin)
      else
        is1beginrange = ns1
        do is = isbegin, isend
          if (is1(is).ne.-1) then
            is1beginrange = min(is1beginrange, is1(is))
          endif
        enddo
      endif      

      return
      end function is1beginrange

      integer function is1endrange(isbegin,isend,ns,ns1,is1)
      implicit none
      integer isbegin, isend, ns, ns1
      integer is1(0:ns-1)
      integer is
      intrinsic max

      if (isbegin.eq.0 .and. isend.eq.ns-1) then
        is1endrange = ns1-1
      else if (isbegin.eq.isend .and. is1(isend).ne.-1) then
        is1endrange = is1(isend)
      else
        is1endrange = -1
        do is = isbegin, isend
          if (is1(is).ne.-1) then
            is1endrange = max(is1endrange, is1(is))
          endif
        enddo
      endif      

      return
      end function is1endrange

      integer function iatmbeginrange(ibegin,iend,natmi,natmi1,iatm1)
      implicit none
      integer ibegin, iend, natmi, natmi1
      integer iatm1(0:natmi-1)
      integer iatm
      intrinsic min

      if (ibegin.eq.0) then
        iatmbeginrange = 0
      else if (ibegin.eq.iend .and. iatm1(ibegin).ne.0) then
        iatmbeginrange = iatm1(ibegin)
      else
        iatmbeginrange = natmi1
        do iatm = ibegin, iend
          if (iatm1(iatm).ne.0) then
            iatmbeginrange = min(iatmbeginrange, iatm1(iatm))
          endif
        enddo
      endif      

      return
      end function iatmbeginrange

      integer function iatmendrange(ibegin,iend,natmi,natmi1,iatm1)
      implicit none
      integer ibegin, iend, natmi, natmi1
      integer iatm1(0:natmi1-1)
      integer iatm
      intrinsic max

      if (iend.eq.0) then
        iatmendrange = 0
      else if (ibegin.eq.iend .and. iatm1(iend).ne.0) then
        iatmendrange = iatm1(ibegin)
      else
        iatmendrange = 0
        do iatm = ibegin, iend
          if (iatm1(iatm).ne.0) then
            iatmendrange = max(iatmendrange, iatm1(iatm))
          endif
        enddo
      endif      

      return
      end function iatmendrange


