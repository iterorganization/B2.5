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

      subroutine b2ytdr (ninp, nout, nx1, ny1)
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
      use b2mod_feedback
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
      use b2mod_b2plot
     , ,only: nx, ny, ns, wklng,
     ,        ixref, iyref,
     ,        jsep, jxa, jxi, nxtl, nxtr,
     ,        nstsi, nlimi, natmi, nmoli, nioni,
     ,        alloc_b2mod_b2plot, dealloc_b2mod_b2plot
      use b2mod_external
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_braeir
      use eirmod_comusr
     , ,only: NSPH, NMASSA, NCHARA
      use eirmod_ccoupl
      use eirmod_eirbra
      use eirmod_extrab25
     , ,only: msrfclsx, mgwtiesx, lkindp, lkindm, lkindi, texta,
     ,        eirene_extrab25_eirpbls_init, eirene_extrab25_wneuinit,
     ,        eirene_extrab25_cleanup
      use eirmod_clogau
#endif
      implicit none
#     include "DIMENSIONS.F"
#ifndef B25_EIRENE
#     include <COUPLE/BRAEIR.F>
#endif

*   ..input arguments (unchanged on exit)
      integer ninp(0:8), nout(0:5), nx1, ny1
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
*     ninp - (0:8) integer array, input.
*     ninp specifies the input unit numbers as follows.
*     ninp(0): formatted; specifies the output mesh size.
*     ninp(1): un*formatted; provides the plasma state.
*     ninp(2): formatted; provides the run switches.
*     ninp(3): un*formatted; provides the (new) geometry.
*     ninp(4): un*formatted; provides the (new) atomic rates.
*     ninp(5): formatted; provides the transport switches.
*     ninp(6): formatted; provides the initial state.
*     ninp(7): formatted; provides the (old) geometry.
*     ninp(8): un*formatted; provides default physics parameters.
*
*     nout - (0:5) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides printed output.
*     nout(1): graphical; provides the transferred plasma state.
*     nout(2): formatted; provides modified run switches.
*     nout(3): formatted; provides modified transport parameters.
*     nout(4): formatted; provides modified initial state.
*     nout(5): formatted; provides the list of new modified files.
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
      integer, allocatable :: idum(:)
      integer ii(1), i1, i2, ik, ikind, inc, nnc
      integer nbc1, ndat, nds(ns-1)
      integer boundary_namelist, neutrals_namelist, use_eirene
      integer numerics_namelist, sources_inputfile, b2stbc_feedback
      integer transport_namelist, transport_inputfile
      integer atomic_physics_rescale_flag
      integer ismain, ismain0, isfeedback, isfb
      integer nsmin, nsmax, ns1min, ns1max
      integer sput_phys_col, sput_phys_model, refl_model
      integer iy_nocoreExB, bc_type16_ref, bc_type20_ref
      integer ixref1, iyref1
      integer ix1, ix2, iy1, iy2
      integer ix11, ix21, iy11, iy21
      integer elm_ix1, elm_ix2
      integer mw_ix1, mw_ix2, mw_jxa, mw_jsep
      integer mw_jxa1, mw_jsep1
      integer jsep1, jxi1, jxa1, nxtl1, nxtr1, iswitch
      integer nncut1, nnreg1(0:2), periodic_bc1
      integer scaling_ix_begin, scaling_ix_end
      integer elm_ix_begin, elm_ix_end
      integer ndim, ndm2, nfitanf, nfitend, ifail
      character lblmn*120, lblgm*120, lblrc*120, lblcp*120, lblps*120,
     .          ligne*120
      character id0*8
      character*8 eirene_format
      character*4 solps_version
      data eirene_format / 'iter' /, solps_version / 'iter' /
      character*10 :: b2fgmtry_oldversion, b2fgmtry_newversion
      data b2fgmtry_oldversion / '00.000.000' /
      character tab
      real (kind=R8) ::
     1 fhe1(-1:nx1,-1:ny1,0:1),
     2 fch1(-1:nx1,-1:ny1,0:1),
     3 fhi1(-1:nx1,-1:ny1,0:1), po1(-1:nx1,-1:ny1),
     5 ti1(-1:nx1,-1:ny1), te1(-1:nx1,-1:ny1), ne1(-1:nx1,-1:ny1),
     6 fch1_32(-1:nx1,-1:ny1,0:1), fch1_52(-1:nx1,-1:ny1,0:1),
     7 fch1_p(-1:nx1,-1:ny1,0:1), w1d(0:ns-1), zpdum(0:ns-1,0:3)
      real (kind=R8), allocatable ::
     1 zamin1(:), zamax1(:), zn1(:), am1(:),
     2 uadia1(:,:,:,:), fna1(:,:,:,:), ua1(:,:,:), na1(:,:,:),
     3 kinrgy1(:,:,:),
     4 sputter_yield_yt(:,:,:,:), sputter_yield2_yt(:,:,:,:),
     . sputter_tmp(:,:), cftmp(:), e_bundle_contr(:,:), epls_map(:,:),
     5 bundle_contr(:,:), bundle_frac(:,:), x1data(:,:,:,:),
     6 s1data(:,:,:,:), t1data(:,:,:,:), eta_hce_mult1(:),
     7 recyc1(:,:), erecyc1(:,:), rcion1(:,:),
     8 phys_sput1(:,:), chem_sput1(:,:),
     9 r(:), f(:), ds(:), dprof(:), felm(:), delm(:),
     * px(:), pf(:), pelm(:),
     * prdsa(:), pra(:), dsa(:), dsa1(:), deta(:),
     * pxdsa(:), pxa(:), xx(:), xf(:), x(:), fx(:), dx(:),
     1 sr1(:), sf1(:), sd1(:), sr2(:), sf2(:), sd2(:),
     2 spr(:), spf1(:), spf2(:), s1d(:,:,:,:),
     1 xr1(:), xf1(:), xd1(:), xr2(:), xf2(:), xd2(:),
     2 xpr(:), xpf1(:), xpf2(:), x1d(:,:,:,:),
     3 l_conn(:), nbini(:),
     4 saved_na_feedback_actuator1(:), na_feedback_target1(:),
     5 na_feedback_alpha1(:), na_feedback_beta1(:),
     6 na_feedback_puff_min1(:), na_feedback_puff_max1(:),
     7 na_feedback_const1(:), na_feedback_overshoot1(:),
     8 na_feedback_time1(:)
      real (kind=R8) :: ds_offset, l_total, l_ref, na0eps
      integer ns1, nsm, nym, natmm, nwallm
      integer i, j, k, ix, iy, is, ib, ireg, len, nsns
      integer ib1, ib2, ifl, js, ks, kmod, kint, itrack
      integer natm1, nspecies1, iss
      integer itrack1, itrack2, itrack3
      integer non_commensurate
      integer is0(0:DEF_NSD-1), is1(0:ns-1)
      integer is1beginrange, is1endrange, iatmbeginrange, iatmendrange
      integer, allocatable :: iatm1(:), ipls1(:), istra1(:)
      integer, allocatable :: istmp(:), ns1data(:,:,:), nx1data(:,:,:)
      integer, allocatable :: mlcmp1(:,:), ifl1(:), ispecies1(:)
      integer, allocatable :: b2espcr1(:), eb2spcr1(:), eb2ibcr(:)
      integer, allocatable :: e_bundle_start(:), e_bundle_end(:)
      integer, allocatable :: n1data(:,:,:), ad1spec(:,:,:), isp1(:)
      integer, allocatable :: lstrascl1(:,:)
      integer, allocatable :: na_feedback_actuator1(:),
     1  na_feedback_choice1(:), na_feedback_option1(:),
     2  na_feedback_ix11(:), na_feedback_ix21(:),
     3  na_feedback_iy11(:), na_feedback_iy21(:),
     4  na_feedback_ib1(:)
      integer bcpos1(nbcd), bcstart1(nbcd), bcend1(nbcd)
      integer bc_list_size1(nbcd)
      integer bc_list_x1(2*(nx1+ny1),nbcd), bc_list_y1(2*(nx1+ny1),nbcd)
      integer rcpos1(nstraid), rcstart1(nstraid), rcend1(nstraid)
      integer rc_list_size1(nstraid)
      integer rc_list_x1(2*(nx1+ny1),nstraid),
     .        rc_list_y1(2*(nx1+ny1),nstraid)
      integer species_start1(nstraid), species_end1(nstraid)
      integer targsp1(nstraid,ntrack)
      integer nxiso3(nnisod), nxiso4(nnisod)
      integer nyiso3(nnisod), nyiso4(nnisod)
      integer ixa, ixb, iya, iyb
      integer nstrai1, nstrat1
#ifdef B25_EIRENE
      integer IUNIN, IUNOT, INDGRD(3), NR1ST, NP2ND, NT3RD
      integer INDPRO(12), NMODE, NAINI
      integer NZADD, NVOLPR, NSURPR, NVOLPL, NSP
      integer NPLINR,NPLOTR,NPLDLR,NPLINP,NPLOTP,NPLDLP,
     .        NPLINT,NPLOTT,NPLDLT
      integer IDIMP, INUMP, IRPTA1, IRPTE1, IRPTA2, IRPTE2, ILEQUI
      integer IRPTA3, IRPTE3
      integer NPHOTI, NPLSI, NNREAC, IPLTA, ILINE
      integer NATMI1, NREAC1, NSTRA1, NPLSI1, NPLSM
      integer IREAD, IEND, IANF, ISTS, IH, I3, L, M, JSS
      integer IATM, IMOL, IION, IPHOT, ICRM, IPLS, IPLSS, ISTRA
      integer MAXRAD, NSIGVI, NSIGSI, NSIGCI, NPPLG, NCRM
      integer NSIGV1, NSIGS1, NSIGC1
      integer NADVI, NCLVI, NALVI, NADSI, NALSI
      integer ISPL, IPROJ, NPLI1, NPLO1, NPLI2, NPLO2
      integer NPLIN2, NPLOT2, NPLDL2
      integer NCHORI, LNTSPEC, LNTSPUT, NLSTADDS, NGWTIES, NSNVI, NPWIN
      integer iside, numpairs
      integer NTL, NFR
      integer IPOS0, IPOS1, IPOS2
      integer NPRCSF, IPRCSF
      integer MAXPOL, MAXTOR, MAXADD
      integer NSPZ1
      integer NSIGI_BGK, NSIGI_COP, NSIGI_SPC
      integer NADMOD, NASMOD, NORMOD, IGEOUSR
      integer NSSIR, NSSIP, IR, IP
      integer, allocatable :: NPOINT(:,:), ILIIN(:), IHELP(:), IGJUM0(:)
      integer, allocatable :: IREAC1(:)
      integer, allocatable :: INDSRC(:)
      integer, allocatable :: NPRTA(:),
     &                        ISRFA(:), ISRTA(:), ISRSA(:), NRCA(:),
     &                        IREACA(:,:), IBULKA(:,:),
     .                        ISCD1A(:,:), ISCD2A(:,:), ISCD3A(:,:),
     .                        ISCD4A(:,:), ISCDEA(:,:)
      integer, allocatable :: NMASSM(:), NCHARM(:), NPRTM(:),
     &                        ISRFM(:), ISRTM(:), ISRSM(:), NRCM(:),
     &                        IREACM(:,:), IBULKM(:,:),
     .                        ISCD1M(:,:), ISCD2M(:,:), ISCD3M(:,:),
     .                        ISCD4M(:,:), ISCDEM(:,:)
      integer, allocatable :: NMASSI(:), NCHARI(:), NPRTI(:), NCHRGI(:),
     &                        ISRFI(:), ISRTI(:), ISRSI(:), NRCI(:),
     &                        IREACI(:,:), IBULKI(:,:),
     .                        ISCD1I(:,:), ISCD2I(:,:), ISCD3I(:,:),
     .                        ISCD4I(:,:), ISCDEI(:,:)
      integer, allocatable :: ISRFPH(:), ISRTPH(:), ISRSPH(:), NRCPH(:),
     &                        IREACPH(:,:), IBULKPH(:,:),
     .                        ISCD1PH(:,:), ISCD2PH(:,:), ISCD3PH(:,:),
     .                        ISCD4PH(:,:), ISCDEPH(:,:)
      integer, allocatable :: NMASSP(:), NCHARP(:), NPRTP(:), NCHRGP(:),
     &                        ISRFP(:), ISRTP(:), ISRSP(:), NRCP(:),
     &                        DENSLIM(:), NRE(:), NHSTP(:), TDMISP(:,:),
     &                        IREACP(:,:), IBULKP(:,:),
     .                        ISCD1P(:,:), ISCD2P(:,:), ISCD3P(:,:),
     .                        ISCD4P(:,:), ISCDEP(:,:)
      integer, allocatable :: NSCRM(:), NIONREAC(:), NRCMBREAC(:),
     &                        TYPCRM(:),
     .                        SCRM(:,:), NNCRM(:,:), PHCRM(:,:)
      integer, allocatable :: ISP(:), NMASSP1(:), NCHARP1(:),
     &                        NPRTP1(:), NCHRGP1(:),
     &                        ISRFP1(:), ISRTP1(:), ISRSP1(:), NRCP1(:),
     &                        DENSLIM1(:), NRE1(:), NHSTP1(:),
     &                        TDMISP1(:,:)
      integer, allocatable :: IPLAA(:), IPLEA(:)
      integer, allocatable :: IPLAS(:), IPLES(:)
      integer, allocatable :: IFLB1(:), LKINDP1(:)
      integer, allocatable :: NPTS(:), NINITL(:), NEMODS(:), NAMODS(:)
      integer, allocatable :: NMINPTS(:), NPTSDEL(:), NRAYEN(:)
      integer, allocatable :: NSPEZ(:)
      integer, allocatable :: NSRFSI(:)
      integer, allocatable :: INDIM(:,:), INSOR(:,:),
     .                        INGRDA(:,:,:), INGRDE(:,:,:)
      integer, allocatable :: NRSOR(:,:), NPSOR(:,:), NTSOR(:,:),
     .                        NBSOR(:,:), NASOR(:,:), NISOR(:,:),
     .                        ISTOR(:,:)
      integer, allocatable :: IGH(:), IIH(:), IGHW(:), IIHW(:),
     .                        IGHC(:,:), IIHC(:,:)
      integer, allocatable :: NFOLA(:), NGENA(:), NHSTA(:)
      integer, allocatable :: NFOLM(:), NGENM(:), NHSTM(:)
      integer, allocatable :: NFOLI(:), NGENI(:), NHSTI(:)
      integer, allocatable :: NFOLPH(:), NGENPH(:), NHSTPH(:)
      integer, allocatable :: IESTMA(:,:), IBGKA(:,:)
      integer, allocatable :: IESTMM(:,:), IBGKM(:,:)
      integer, allocatable :: IESTMI(:,:), IBGKI(:,:)
      integer, allocatable :: IESTMPH(:,:), IBGKPH(:,:)
      integer, allocatable :: IVLSF(:),
     .                        ISCLS(:), ISCLT(:), ISCL1(:),
     .                        ISCL2(:), ISCL3(:), ISCLB(:),
     .                        ISCLA(:)
      real (kind=R8) :: rdum(0:17), RMN, RMX, ALLOC, MPTS
      real (kind=R8) :: zmin1, zmin2, zmax1, zmax2
      real (kind=R8), allocatable :: XPOL(:,:), YPOL(:,:)
      real (kind=R8), allocatable :: VOLCOR(:), VOLADD(:)
      real (kind=R8), allocatable :: RLB(:)
      real (kind=R8), allocatable :: EELECA(:,:), EBULKA(:,:),
     .           ESCD1A(:,:), ESCD2A(:,:), FREACA(:,:), FLDLMA(:,:)
      real (kind=R8), allocatable :: EELECM(:,:), EBULKM(:,:),
     .           ESCD1M(:,:), ESCD2M(:,:), FREACM(:,:)
      real (kind=R8), allocatable :: EELECI(:,:), EBULKI(:,:),
     .           ESCD1I(:,:), ESCD2I(:,:), FREACI(:,:)
      real (kind=R8), allocatable :: EELECPH(:,:),EBULKPH(:,:),
     .           ESCD1PH(:,:),ESCD2PH(:,:),FREACPH(:,:),FLDLMPH(:,:)
      real (kind=R8), allocatable :: EELECP(:,:), EBULKP(:,:),
     .           ESCD1P(:,:), ESCD2P(:,:), FREACP(:,:)
      real (kind=R8), allocatable ::
     .           TDMBOLTZ(:), TDMDELTAE(:), TDM_CORONA(:)
      real (kind=R8), allocatable ::
     .           TDMBOLTZ1(:), TDMDELTAE1(:), TDM_CORONA1(:)
      real (kind=R8), allocatable ::
     .  TIO(:),TIl(:),TI2(:),TI3(:),TI4(:),TI5(:)
      real (kind=R8), allocatable ::
     .  DI0(:),DI1(:),DI2(:),DI3(:),DI4(:),DI5(:)
      real (kind=R8), allocatable ::
     .  VX0(:),VX1(:),VX2(:),VX3(:),VX4(:),VX5(:)
      real (kind=R8), allocatable ::
     .  VY0(:),VY1(:),VY2(:),VY3(:),VY4(:),VY5(:)
      real (kind=R8), allocatable ::
     .  VZ0(:),VZ1(:),VZ2(:),VZ3(:),VZ4(:),VZ5(:)
      real (kind=R8), allocatable ::
     .  P0(:),P1(:),P2(:),P3(:),P4(:),P5(:)
      real (kind=R8), allocatable :: DATD(:), DATD1(:)
      real (kind=R8), allocatable :: DPLD(:), DPLD1(:)
      real (kind=R8), allocatable :: FCTE1(:)
      real (kind=R8), allocatable :: BMASS1(:)
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
      real (kind=R8), allocatable :: RAYFRAC(:)
      real (kind=R8), allocatable :: SCALV(:)
      real (kind=R8), allocatable :: XTRIAN(:), YTRIAN(:)
#endif
      real (kind=R8) ::
     *  neutral_rescale, BoRiS, divheat, rswitch,
     .  scaling_width, scaling_strength, ixfb, iyfb, ixfb1, iyfb1,
     .  elm_time_begin, elm_time_end, elm_time_period,
     .  sources_time_mod, sources_time_switch,
     1  zamin0(0:DEF_NSD-1), zamax0(0:DEF_NSD-1),
     2  zn0(0:DEF_NSD-1), am0(0:DEF_NSD-1),
     4  rescale_ra(0:nsdecl-1), rescale_sa(0:nsdecl-1),
     5  rescale_qa(0:nsdecl-1), rescale_cx(0:nsdecl-1),
     6  rescale_rd(0:nsdecl-1), rescale_br(0:nsdecl-1)
      real (kind=R8) :: recyceir1(nstraid), volrecstart1(nstraid),
     &  userfluxparm1(nstraid,2)
      character cnamip*80, cvalip*80
      character*256 sources_filename, new_filename, filename
      character*1 crcstra1(nstraid)
      character*36 hlp_frm, hlp_frp, hlp_frb
      character*42 fmt0
      parameter (fmt0='(0p,1x,'''''''',i3,'''''''',1x,4(f9.3,1x))')
      character*500 line
#ifdef B25_EIRENE
      character*80 ZEILE, ZEILEN(100)
      character*80, allocatable :: REAC_LINE(:)
      character*24 TXTSPC, TXTUNT, TEXTLA
      character*72 TXTMOL, TXTION, TXTPLS, TXTPHOT, TXTCRM
      character*72 TITLE5, CHRTAL
      character*74 PATH, FILE, RFILNM
      character HANDLE*6
      character FILNAM*8,H123*4,REAC*9,REAC2*50,CRC*3
      character ELNAME*2,REACSTR*11
      character*8 varname
      character*8, allocatable ::
     .  TEXTM(:), TEXTI(:), TEXTP(:), TEXTPH(:)
      character*8, allocatable :: TEXTP1(:)
      character*72, allocatable :: TXTSOU(:)
      character*10, allocatable :: CDENMODEL(:), CDENMODEL1(:)
      character*6, allocatable :: TDMFNAME(:,:), TDMFNAME1(:,:)
      character*4, allocatable :: TDMH2(:,:), TDMH21(:,:)
      character*9, allocatable :: TDMREACTION(:,:), TDMREACTION1(:,:)
      character*3, allocatable :: TDMCR(:,:), TDMCR1(:,:)
#endif
      character*10, allocatable :: species_name(:)
      logical ldum(0:60), found, ex, l_chem_sput, is_core_boundary
      logical chemsp1(nstraid), savestratum(0:nstraid)
      logical, allocatable :: is_neutral1(:)
      logical, allocatable :: savesigv(:), savesigs(:), savesigc(:)
#ifdef B25_EIRENE
      logical, allocatable :: NLAVRP(:), NLAVRT(:), NLSYMP(:), NLSYMT(:)
      logical, allocatable :: NLRAY(:)
      logical, allocatable :: NLATM(:), NLMOL(:), NLION(:), NLPLS(:),
     .                        NLPHOT(:)
      logical, allocatable :: NLPNT(:), NLLNE(:), NLSRF(:),
     .                        NLVOL(:), NLCNS(:),
     .                        NLPRCA(:), NLPRCM(:), NLPRCI(:),
     .                        NLPRCPH(:), NLPRCA1(:)
      logical, allocatable :: TRCSRC(:), TRCSR1(:)
      logical, allocatable :: nlpbls(:), nlpbl1(:)
      logical, allocatable :: done(:)
      logical LRPSCUT
      logical in_bundle, charges_match, skip_reac
#endif
      logical todo, same_bundling, first_pass
      logical wall_namelist_used
      logical region_flags(NNREGMAX,nkind_coeff), no_pflux
      logical poloidal_scaling, prdone, axdone
      logical clash(ns-1), shared(0:ns-1)
      logical track_chem_sput1(ntrack)
*   ..procedures
      integer lnblnk, convert_position, len_of_digits, transform_tally
      integer conformal_convert, convert_eirene, b2idez
      real (kind=R8) :: cr, cz, cr1, cz1, rratio
      logical open_file, streql, is_comment
      external subini, subend, xertst, xerrab, cfruch, cfrure, cfwuin,
     &  cfwuch, cfwure, ipgeti, ipgetc, lnblnk, convert_position,
     &  len_of_digits, open_file, get_jsep, conformal_convert,
     &  convert_eirene, b2idez, streql, get_nxt, b2yt_nxt, b2yt_ngread,
     &  is_comment
      external b2rups, b2rurc, b2wups, b2xtps, b2rflb, b2ytgm, rratio
      external b2ahx0, b2aix0, b2rfcp, b2wfcp, b2xvcp, b2rugm, cfverr
      external E01BEF, E01BFF
      intrinsic nint, mod, min0, max, index, iabs, sqrt
*   ..initialisation
      data non_commensurate/0/
      data boundary_namelist/0/, neutrals_namelist/0/, use_eirene/0/
      data numerics_namelist/0/, sources_inputfile/0/
      data transport_namelist/0/, transport_inputfile/0/
      data atomic_physics_rescale_flag/0/
      data itrack1/-1/, itrack2/-1/, itrack3/-1/
      data ismain/1/, ismain0/0/, isfeedback/0/, b2stbc_feedback/0/
      data sput_phys_col/3/, sput_phys_model/1/, refl_model/1/
      data iy_nocoreExB/-2/, bc_type16_ref/1/, bc_type20_ref/1/
      data mw_jxa/-1/, mw_jsep/-1/
      data ix1/-1/, ix2/-1/, iy1/-1/, iy2/-1/
      data elm_ix1/-1/, elm_ix2/-1/, mw_ix1/-1/, mw_ix2/-1/
      data BoRiS/0.0_R8/, na0eps/0.0_R8/
*  ..namelists
      namelist /profile/ nsdata, nxdata, sdata, xdata, divheat,
     . sources_time_mod, sources_time_switch, sources_filename
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
*     The sign of the parallel velocity indicates whether it is flowing
*     in the same direction as the local poloidal magnetic field.
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
*     (ix,iy) cell and its left neighbor.
*     For ix.eq.-1, fnax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fnay(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix,iy) cell and its bottom neighbor.
*     For iy.eq.-1, fnay(-1:nx,iy,0:ns-1) will hold 0.
*
*     fhe - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhex(,)=fhe(,,0) and fhey(,)=fhe(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhex(ix,iy) specifies the electron
*     heat flux through the face between the (ix,iy) cell and its
*     left neighbor. For ix.eq.-1, fhex(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhey(ix,iy) specifies the electron
*     heat flux through the face between the (ix,iy) cell and its bottom
*     neighbor. For iy.eq.-1, fhey(-1:nx,iy) will hold 0.
*
*     fhi - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhix(,)=fhi(,,0) and fhiy(,)=fhi(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhix(ix,iy) specifies the all atom
*     heat flux through the face between the (ix,iy) cell and its
*     left neighbor. For ix.eq.-1, fhix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhiy(ix,iy) specifies the all atom
*     heat flux through the face between the (ix,iy) cell and its
*     bottom neighbor. For iy.eq.-1, fhiy(-1:nx,iy) will hold 0.
*
*     fch - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchx(,)=fch(,,0) and fchy(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchx(ix,iy) specifies the electric
*     current through the face between the (ix,iy) cell and its left
*     neighbor. For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchy(ix,iy) specifies the electric
*     current through the face between the (ix,iy) cell and its bottom
*     neighbor. For iy.eq.-1, fchy(-1:nx,iy) will hold 0.
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
     &             1.le.ninp(3).and.1.le.ninp(4).and.1.le.ninp(5).and.
     &             1.le.ninp(6).and.1.le.ninp(7).and.1.le.ninp(8),
     &             'faulty argument ninp')
      call xertst (1.le.nout(0).and.1.le.nout(1).and.1.le.nout(2).and.
     &             1.le.nout(3).and.1.le.nout(4).and.1.le.nout(5),
     &             'faulty argument nout')
*   ..test dimensions
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')
      call xertst (0.le.nx1.and.0.le.ny1, 'faulty argument nx1, ny1')
      write(*,*) 'nx, nx1 = ', nx, nx1
      write(*,*) 'ny, ny1 = ', ny, ny1
#ifndef B25_EIRENE
      allocate(idum(0:12))
#else
      allocate(idum(0:max(12,msrfclsx,nlimps,mgwtiesx)))
#endif
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
      call b2yt_nxt (nx1, nxtl1, nxtr1)
*   ..read and echo code internal parameters from b2yt.dat
      write (nout(0),'(/2x,a)')
     . 'non-default internal parameters from b2yt.dat:'
    3 continue
       read (ninp(0),'(a)',end=4,err=91) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         write (nout(0),'(4x,a,2x,a)') cnamip, cvalip
         call ipsetc (cnamip, cvalip)
       endif
       goto 3
    4  continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'
#ifdef B25_EIRENE
      call ipgetc ('eirene_format', eirene_format)
      write(*,*) 'Eirene input.dat format entered : ', eirene_format
      if (streql(eirene_format,'iter')) then
        write(*,*)
     .   'Assuming standard SOLPS-ITER input.dat format (default)'
      else if (streql(eirene_format,'old')) then
        write(*,*) 'Assuming "old" Eirene_96 input.dat format'
      else if (streql(eirene_format,'new')) then
        write(*,*) 'Assuming "new" Eirene_99 input.dat format'
      else if (streql(eirene_format,'facelift')) then
        write(*,*) 'Assuming Eirene_facelift input.dat format'
      else if (streql(eirene_format,'juelich')) then
        write(*,*) 'Assuming Juelich Eirene input.dat format'
      else
        write(*,*) 'Accepted values are (case-insensitive): '
        write(*,*) '"ITER", "old", "new", "facelift", "Juelich"'
        call xerrab ('Unrecognized Eirene input.dat format!')
      endif
      call ipgetc ('solps_version', solps_version)
      write(*,*) 'Incoming SOLPS version entered : ', solps_version
      if (streql(solps_version,'iter')) then
        write(*,*)
     .   'Assuming standard SOLPS-ITER fort.44 format (default)'
      else if (streql(solps_version,'4.3')) then
        write(*,*) 'Assuming SOLPS4.3 fort.44 format'
      else if (streql(solps_version,'5.0')) then
        write(*,*) 'Assuming SOLPS5.0 fort.44 format'
      else if (streql(solps_version,'5.1')) then
        write(*,*) 'Assuming SOLPS5.1 fort.44 format'
      else if (streql(solps_version,'5.2')) then
        write(*,*) 'Assuming SOLPS5.2 fort.44 format'
      else
        write(*,*) 'Accepted values are (case-insensitive): '
        write(*,*) '"ITER", "4.3", "5.0", "5.1", "5.2"'
        call xerrab ('Unrecognized SOLPS version format!')
      endif
*   ..import parameters from Eirene input.dat file
      inquire(file='input.dat',exist=ex)
      if (ex) then
        open(unit=45,file='input.dat')
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
        write(*,*) '  NPLS = ',NPLS
        write(*,*) 'NLIMPS = ',NLIMPS
      else
        ndx=DEF_NXD
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
        N2ND=DEF_NXD+1+DEF_NCUT+DEF_ISOEXTRA
        N3RD=DEF_NPHID
        CALL EIRENE_SET_PARMMOD(1)
        CALL EIRENE_ALLOC_CLOGAU
      endif
      NATMI = NATM
      NPLSI = NPLS
      NMOLI = NMOL
      NIONI = NION
      CALL EIRENE_ALLOC_BRAEIR(NDX,NDY,NFL,ifoff)
      CALL EIRENE_ALLOC_EIRBRA(NDX,NDY,NFL,NSTRA,ifoff)
      CALL EIRENE_ALLOC_CCOUPL(1)
      CALL EIRENE_SETCON
      CALL EIRENE_EXTRAB25_WNEUINIT
      CALL EIRENE_extrab25_eirpbls_init(1,natm,nmol,nion,npls,nstra,nfl)
#else
      ex = .false.
#endif

* ..computation
*   ..read plasma state file label
      call cfruch (ninp(1), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(1), b2fstate_version,
     &  ns, zamin, zamax, zn, am, .true.)
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
      if (non_commensurate.eq.1 .or. nx.ne.nx1 .or. ny.ne.ny1) then
*   ..check for and open old geometry file
        call xertst (open_file(70,'b2fgmtry.old',''),
     .   'Geometry translation requires b2fgmtry.old file')
*   ..obtain old geometry version number
        call cfverr (70, b2fgmtry_oldversion)
        if (b2fgmtry_oldversion.ne.b2fgmtry_version) then
          b2fgmtry_newversion=b2fgmtry_version
          b2fgmtry_version=b2fgmtry_oldversion
        endif
*   ..obtain nx, ny from old geometry file
        call cfruin (70, 2, idum, 'nx,ny')
        call xertst (nx.eq.idum(0),
     &  'inconsistent nx definition between b2fgmtry.old and b2fstate!')
        call xertst (ny.eq.idum(1),
     &  'inconsistent ny definition between b2fgmtry.old and b2fstate!')
*   ..read header for old geometry
        call cfruch (70, 120, lblgm, 'label')
        call cfruin (70, 1, idum, 'isymm')
*   ..read old geometry and metric
        call b2rugm (70, nx, ny, crx, cry, fpsi,
     &   ffbz, bb, vol, hx, hy, qz, qc, gs, pbs, wbbl, wbbr, wbbv)
        if (non_commensurate.eq.1) then
        call xertst (nnreg(0).eq.nnreg1(0) .and. nncut.eq.nncut1 .and.
     &   periodic_bc.eq.periodic_bc1,
     &   'Non-commensurate translation requires the same topology!')
        if (nncut.eq.2) call xertst(
     .   (topcut(1).eq.topcut(2) .and. topcut1(1).eq.topcut1(2)) .or.
     .   (topcut(1).lt.topcut(2) .and. topcut1(1).lt.topcut1(2)) .or.
     .   (topcut(1).gt.topcut(2) .and. topcut1(1).gt.topcut1(2)),
     .   'Non-commensurate translation requires the same DN topology!')
        end if
        close(70)
        if (b2fgmtry_oldversion.ne.b2fgmtry_newversion) then
          b2fgmtry_version=b2fgmtry_oldversion
        endif
      else if (nx.eq.nx1 .and. ny.eq.ny1) then
*   ..easy case where we keep the same grid geometry
*   ..simply re-read the file to populate the b2mod_geo arrays
        rewind(ninp(3))
        call cfverr (ninp(3), b2fgmtry_oldversion)
        call cfruin (ninp(3), 2, idum, 'nx,ny')
        call cfruch (ninp(3), 120, lblgm, 'label')
        call cfruin (ninp(3), 1, idum, 'isymm')
        call b2rugm (ninp(3), nx, ny, crx, cry, fpsi,
     &   ffbz, bb, vol, hx, hy, qz, qc, gs, pbs, wbbl, wbbr, wbbv)
      endif
      call get_jsep(nx,ny,jxi,jxa,jsep)
      call get_nxt(nx,nxtl,nxtr)
      call alloc_b2mod_b2plot(nx,ny,ns,ns_ext,wklng,
     1 natm,nmol,nion,nlimps)
      nsm=max(ns,ns1)
      nym=max(ny,ny1)
      call alloc_input_profile(nym,nsm)
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
      allocate(bundle_contr(0:ns-1,0:DEF_NSD-1))
      allocate(b2espcr1(0:ns1-1))
      allocate(eb2ibcr(npls))
      allocate(ispecies1(natm))
      allocate(ifl1(ns1))
      allocate(cftmp(0:nsm-1))
      allocate(istmp(0:nsm-1))
      allocate(iatm1(natm))
      allocate(ipls1(npls))
      bundle_frac=0.0_R8
      bundle_contr=0.0_R8
      ispecies1=-1
      nsmin=0
      nsmax=ns
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
        write(hlp_frm,'(a,i3,a)') '(3x,',ns,'(1x,i4))'
        write(*,hlp_frm) (is, is = 0, ns-1)
        write(hlp_frm,'(a,i3,a)') '(i3,',ns,'f5.2)'
        do iss = 0, ns1-1
          write(*,hlp_frm) iss, (bundle_frac(iss,is), is = 0, ns-1)
        end do
        nspecies1 = 1
        ib = 1
        do while (ib.lt.ns1)
          if (nint(zn1(ib)).ne.nint(zn1(ib-1)) .or.
     &        nint(am1(ib)).ne.nint(am1(ib-1))) then
            nspecies1=nspecies1+1
          end if
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
      natm1 = nspecies1
      allocate(eb2spcr1(natm1))
      allocate(isp1(natm1))
      allocate(mlcmp1(natm1,nmol))
      eb2spcr1=-1
      b2espcr1=-1
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
* ..b2mn.dat, b2.*.parameters, input.dat, b2ah.dat, b2ai.dat

*   ..read default physics common
      call cfruch (ninp(8), 120, lblcp, 'label')
      call b2rucp (ninp(8), b2fparam_version, nx, ns)
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
       read (ninp(2),'(a)',end=2,err=91) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         write (nout(0),'(4x,a,2x,a)') cnamip, cvalip
         call ipsetc (cnamip, cvalip)
       endif
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
      call ipgeti ('b2stbc_feedback', b2stbc_feedback)
      call ipgeti ('b2mndr_atomic_physics_rescale',
     .              atomic_physics_rescale_flag)
      call ipgeti ('b2mndr_eirene', use_eirene)
      call ipgetr ('b2mndr_na0eps', na0eps)
      call ipgeti ('b2mndr_ismain', ismain)
      call ipgeti ('b2mwti_ismain0', ismain0)
      call ipgeti ('b2news_nsmin', nsmin)
      call ipgeti ('b2news_nsmax', nsmax)
      ixfb=rratio(5*nx,8)
      iyfb=rratio(ny,2)
      isfb=ismain
      call ipgetr ('b2mndr_ixfb', ixfb)
      call ipgetr ('b2mndr_iyfb', iyfb)
      call ipgeti ('b2mndr_isfb', isfb)
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
      inquire(file='../b2.boundary.parameters',exist=todo)
      if (boundary_namelist.eq.1.or.todo)
     > call alloc_b2mod_boundary(nx, ny, ns, nbcd, nnisod)
      inquire(file='../b2.sputter_save.parameters',
     1 exist=sputter_namelist_used)
      inquire(file='../b2.wall_save.parameters',
     1 exist=wall_namelist_used)
      if (sputter_namelist_used.or.wall_namelist_used)
     > call alloc_b2mod_sputter(nx,ny,ns)
      if (transport_namelist.ge.1) call alloc_transport_namelist(ns)
      if (numerics_namelist.ge.1) call alloc_b2mod_numerics(nx, ny, ns)
      inquire(file='../b2.neutrals.parameters',exist=todo)
      inquire(file='../b2.feedback_save.parameters',exist=found)
      todo = todo .or. found
      inquire(file='../b2.feedback_control.parameters',exist=found)
      todo = todo .or. found
      if (neutrals_namelist.eq.1.or.use_eirene.ne.0.or.todo) then
        call alloc_b2mod_neutrals
     .   (nx, ny, ns, natm, nstraid, nsts, sput_dst)
        call read_neutrals_namelist(nx, ny, ns, .false.)
        call xertst (nstrat.le.nstraid,
     .   'Increase DEF_NSTRA in DIMENSIONS.F')
        call alloc_b2mod_diag(nx1,ny1,ns,nnreg(0),nstrat,natm)
      endif
*   ..modify b2mn.dat parameters
      ismain = is1(ismain)
      ismain0 = is1(ismain0)
      isfeedback = is1(isfeedback)
      ixfb1 = ixfb*rratio(nx1,nx)
      iyfb1 = iyfb*rratio(ny1,ny)
      isfb = is1(isfb)
      if (nsmin.gt.ns) then
        ns1min = ns1
      else if (nsmin.lt.0) then
        ns1min = 0
      else if (nsmin.eq.0 .and. nsmax.eq.ns) then
        ns1min = 0
      else if (is1(nsmin).ne.-1) then
        ns1min = is1(nsmin)
      else
        is = nsmin
        do while (is1(is).eq.-1 .and. is.lt.ns)
          is = is+1
        end do
        if (is.lt.ns) then
          ns1min = is1(is)
        else
          ns1min = ns1
        endif
      endif
      if (nsmax.gt.ns) then
        ns1max = ns1
      else if (nsmax.lt.0) then
        ns1max = 0
      else if (nsmin.eq.0 .and. nsmax.eq.ns) then
        ns1max = ns1
      else if (is1(nsmax).ne.-1) then
        ns1max = is1(nsmax)
      else
        is = nsmax
        do while (is1(is).eq.-1 .and. is.ge.0)
          is = is-1
        end do
        if (is.ge.0) then
          ns1max = is1(is)
        else
          ns1max = -1
        endif
      endif
      if (ns1max.lt.ns1min) then
        is = ns1min
        ns1min = ns1max
        ns1max = is
        write(6,*) 'Interval [nsmin:nsmax] flipped with new ordering!'
      end if
      if (nsmin.ne.0 .or. nsmax.ne.ns) then
        do is = max(nsmin,0), min(nsmax,ns)-1
          found = is1(is).ge.ns1min .and. is1(is).le.ns1max-1
          if (.not.found) write(6,*) 'Caution: the interval '//
     .     '[nsmin:nsmax-1] does not map continuously!'
        end do
        do is = ns1min, ns1max-1
          js = 0
          do while (js.lt.ns .and. is1(js).ne.is)
            js = js+1
          end do
          if (js.eq.ns .or. js.lt.nsmin .or. js.ge.nsmax) then
            write(6,'(a,i3,a)') 'Species ',is,
     .       ' from new ordering does not '//
     .       'map into old interval [nsmin:nsmax-1] but is'
            write(6,*) 'included in new interval [ns1min:ns1max-1]!'
            write(6,*) 'Is that what you want?'
          end if
        end do
      end if
      itrack = 0
      l_chem_sput = .false.
      if (sput_dst.ge.0) then
        if (is1(sput_dst).ge.0) then
          sput_dst = is1(sput_dst)
          itrack1 = 1
          itrack = 1
          l_chem_sput = nint(zn(sput_dst)).eq.6
        else
          sput_dst = -1
          itrack1 = 0
        end if
      end if
      if (sput_dst2.ge.0) then
        if (is1(sput_dst2).ge.0) then
          if (sput_dst.lt.0) then
            sput_dst = is1(sput_dst2)
            sput_dst2 = -1
            itrack2 = 1
            itrack = 1
          else
            sput_dst2 = is1(sput_dst2)
            itrack2 = 2
            itrack = 2
          endif
        else
          sput_dst2 = -1
          itrack2 = 0
        endif
      endif
      if (sput_dst3.ge.0) then
        if (is1(sput_dst3).ge.0) then
          if (sput_dst.lt.0 .and. sput_dst2.lt.0) then
            sput_dst = is1(sput_dst3)
            itrack3 = 1
            sput_dst2 = -1
            sput_dst3 = -1
            itrack = 1
          else if (sput_dst2.lt.0) then
            sput_dst2 = is1(sput_dst3)
            itrack3 = 2
            sput_dst3 = -1
            itrack = 2
          else
            sput_dst3 = is1(sput_dst3)
            itrack3 = 3
            itrack = 3
          endif
        else
          sput_dst3 = -1
          itrack3 = 0
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
      read (ninp(2),'(a)',end=6,err=91) line
      if (.not.is_comment(line)) then
        read (line,*) cnamip, cvalip
      else if (.not.streql(line,' ')) then
        i = index(ligne,'''')
        if (i.gt.0) read (line(i:132),*) cnamip, cvalip
      endif
*   ..integer switches
      todo = .true.
      ISWITCHES: select case (cnamip)
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
        case('b2news_nsmin')
          iswitch = ns1min
        case('b2news_nsmax')
          iswitch = ns1max
        case('b2mndr_isfb')
          iswitch = isfb
        case('b2trcl_conductive limit')
          cnamip = 'b2trcl_conductive_limit'
          todo = .false.
        case default ISWITCHES
          todo = .false.
      end select ISWITCHES
      if (todo) then
        write(hlp_frm,'(a,i1,a)') '(i',len_of_digits(iswitch),')'
        write(cvalip,hlp_frm) iswitch
      end if
*   ..real switches
      todo = .true.
      RSWITCHES: select case (cnamip)
        case('b2mndr_ixfb')
          rswitch = ixfb1
        case('b2mndr_iyfb')
          rswitch = iyfb1
        case('b2mndr_na0eps')
          if (na0eps.lt.1.0e8_R8) then
            cnamip = 'b2mndr_na_min'
          else
            cnamip = 'b2mndr_na_new'
          end if
          todo = .false.
        case default RSWITCHES
          todo = .false.
      end select RSWITCHES
      if (todo) write(cvalip,'(1p,e10.3)') rswitch
      if (.not.is_comment(line)) then
        write(nout(2),'(a1,a,a,a,a1)')
     .   '''',trim(cnamip),''''//tab//'''',trim(cvalip),''''
      else if (.not.streql(line,' ')) then
        write(nout(2),'(a,a2,a,a,a,a1)') line(1:1),
     .   ' ''',trim(cnamip),''''//tab//'''',trim(cvalip),''''
      else
        write(nout(2),1028)
      endif
      goto 5
    6 continue

*  ..rewrite boundary namelist
      inquire(file='../b2.boundary.parameters',exist=todo)
      new_filename = '../b2.boundary.parameters.2'
      inquire(file=trim(new_filename),exist=found)
      if (found) write(6,*) trim(new_filename), ' already present!'
      do while((boundary_namelist.eq.1.or.todo).and..not.found)
        call read_b2mod_boundary_namelist(nx, ny, ns, .false.)
        call xertst (nbc.le.nbcd, 'increase size in b2ytdr for nbcd')
        call xertst (nniso.le.nnisod,
     .   'increase size in b2ytdr for nnisod')
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
            if (bcpos1(ib).eq.-2) then
              bc_list_size1(ib) = 0
              do is = 1, bc_list_size(ib)
                ix1 = conformal_convert(non_commensurate,
     .           1,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,1)
                ix2 = conformal_convert(non_commensurate,
     .           1,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,3)
                iy1 = conformal_convert(non_commensurate,
     .           2,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,2)
                do i = ix1, ix2
                  j = 0
                  found = .false.
                  do while (j.le.bc_list_size1(ib) .and. .not.found)
                    j = j + 1
                    found = i.eq.bc_list_x1(j,ib) .and.
     &                    iy1.eq.bc_list_y1(j,ib)
                  end do
                  if (.not.found) then
                    bc_list_size1(ib) = bc_list_size1(ib) + 1
                    bc_list_x1(bc_list_size1(ib),ib) = i
                    bc_list_y1(bc_list_size1(ib),ib) = iy1
                  end if
                end do
              end do
            end if
          else if (bcchar(ib).eq.'N') then
            bcpos1(ib) = conformal_convert(non_commensurate,
     .       2,bcstart(ib),bcpos(ib),nx,ny,nx1,ny1,0)
            bcstart1(ib) = conformal_convert(non_commensurate,
     .       1,bcstart(ib),bcpos(ib),nx,ny,nx1,ny1,1)
            bcend1(ib) = conformal_convert(non_commensurate,
     .       1,bcend(ib),bcpos(ib),nx,ny,nx1,ny1,3)
            if (bcpos1(ib).eq.-2) then
              bc_list_size1(ib) = 0
              do is = 1, bc_list_size(ib)
                ix1 = conformal_convert(non_commensurate,
     .           1,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,1)
                ix2 = conformal_convert(non_commensurate,
     .           1,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,3)
                iy1 = conformal_convert(non_commensurate,
     .           2,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,0)
                do i = ix1, ix2
                  j = 0
                  found = .false.
                  do while (j.le.bc_list_size1(ib) .and. .not.found)
                    j = j + 1
                    found = i.eq.bc_list_x1(j,ib) .and.
     &                    iy1.eq.bc_list_y1(j,ib)
                  end do
                  if (.not.found) then
                    bc_list_size1(ib) = bc_list_size1(ib) + 1
                    bc_list_x1(bc_list_size1(ib),ib) = i
                    bc_list_y1(bc_list_size1(ib),ib) = iy1
                  end if
                end do
              end do
            end if
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
            if (ib2.gt.ib1 .and. bcpos1(ib).ne.-2) then  ! must duplicate the X b.c.
              nbc1 = nbc1 + ib2-ib1
              call xertst (nbc1.le.nbcd,
     .         'increase size in b2ytdr for nbcd')
              do is = nbc1-(ib2-ib1), nbc1
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
                lcbs(is) = 0
              end do
            else if (bcpos1(ib).eq.-2) then
              bc_list_size1(ib) = 0
              do is = 1, bc_list_size(ib)
                ix1 = conformal_convert(non_commensurate,
     .           1,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,1)
                ix2 = conformal_convert(non_commensurate,
     .           1,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,3)
                iy1 = conformal_convert(non_commensurate,
     .           2,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,0)
                iy2 = conformal_convert(non_commensurate,
     .           2,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,2)
                do i = ix1, ix2
                  do j = iy1, iy2
                    k = 0
                    found = .false.
                    do while (k.le.bc_list_size1(ib) .and. .not.found)
                      k = k + 1
                      found = i.eq.bc_list_x1(k,ib) .and.
     &                        j.eq.bc_list_y1(k,ib)
                    end do
                    if (.not.found) then
                      bc_list_size1(ib) = bc_list_size1(ib) + 1
                      bc_list_x1(bc_list_size1(ib),ib) = i
                      bc_list_y1(bc_list_size1(ib),ib) = j
                    end if
                  end do
                end do
              end do
            end if
          else if (bcchar(ib).eq.'W') then
            bcpos1(ib) = conformal_convert(non_commensurate,
     .       1,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,3)
            bcstart1(ib) = conformal_convert(non_commensurate,
     .       2,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,0)
            bcend1(ib) = conformal_convert(non_commensurate,
     .       2,bcpos(ib),bcend(ib),nx,ny,nx1,ny1,2)
            if (bcpos1(ib).eq.-2) then
              bc_list_size1(ib) = 0
              do is = 1, bc_list_size(ib)
                ix1 = conformal_convert(non_commensurate,
     .           1,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,3)
                iy1 = conformal_convert(non_commensurate,
     .           2,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,0)
                iy2 = conformal_convert(non_commensurate,
     .           2,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,2)
                do j = iy1, iy2
                  i = 0
                  found = .false.
                  do while (i.le.bc_list_size1(ib) .and. .not.found)
                    i = i + 1
                    found = j.eq.bc_list_y1(j,ib) .and.
     &                    ix1.eq.bc_list_x1(j,ib)
                  end do
                  if (.not.found) then
                    bc_list_size1(ib) = bc_list_size1(ib) + 1
                    bc_list_x1(bc_list_size1(ib),ib) = ix1
                    bc_list_y1(bc_list_size1(ib),ib) = j
                  end if
                end do
              end do
            end if
          else if (bcchar(ib).eq.'E') then
            bcpos1(ib) = conformal_convert(non_commensurate,
     .       1,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,1)
            bcstart1(ib) = conformal_convert(non_commensurate,
     .       2,bcpos(ib),bcstart(ib),nx,ny,nx1,ny1,0)
            bcend1(ib) = conformal_convert(non_commensurate,
     .       2,bcpos(ib),bcend(ib),nx,ny,nx1,ny1,2)
            if (bcpos1(ib).eq.-2) then
              bc_list_size1(ib) = 0
              do is = 1, bc_list_size(ib)
                ix1 = conformal_convert(non_commensurate,
     .           1,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,1)
                iy1 = conformal_convert(non_commensurate,
     .           2,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,0)
                iy2 = conformal_convert(non_commensurate,
     .           2,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,2)
                do j = iy1, iy2
                  i = 0
                  found = .false.
                  do while (i.le.bc_list_size1(ib) .and. .not.found)
                    i = i + 1
                    found = j.eq.bc_list_y1(j,ib) .and.
     &                    ix1.eq.bc_list_x1(j,ib)
                  end do
                  if (.not.found) then
                    bc_list_size1(ib) = bc_list_size1(ib) + 1
                    bc_list_x1(bc_list_size1(ib),ib) = ix1
                    bc_list_y1(bc_list_size1(ib),ib) = j
                  end if
                end do
              end do
            end if
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
              call xertst (nbc1.le.nbcd,
     .         'increase size in b2ytdr for nbcd')
              do is = nbc1-(ib2-ib1), nbc1
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
                lcbs(is) = 0
              end do
            else if (bcpos1(ib).eq.-2) then
              bc_list_size1(ib) = 0
              do is = 1, bc_list_size(ib)
                ix1 = conformal_convert(non_commensurate,
     .           1,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,1)
                ix2 = conformal_convert(non_commensurate,
     .           1,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,3)
                iy1 = conformal_convert(non_commensurate,
     .           2,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,0)
                iy2 = conformal_convert(non_commensurate,
     .           2,bc_list_x(is,ib),bc_list_y(is,ib),nx,ny,nx1,ny1,2)
                do i = ix1, ix2
                  do j = iy1, iy2
                    k = 0
                    found = .false.
                    do while (k.le.bc_list_size1(ib) .and. .not.found)
                      k = k + 1
                      found = i.eq.bc_list_x1(k,ib) .and.
     &                        j.eq.bc_list_y1(k,ib)
                    end do
                    if (.not.found) then
                      bc_list_size1(ib) = bc_list_size1(ib) + 1
                      bc_list_x1(bc_list_size1(ib),ib) = i
                      bc_list_y1(bc_list_size1(ib),ib) = j
                    end if
                  end do
                end do
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
        do ib = 1, nbc1
          is_core_boundary = .false.
          ib1 = 1
          do while (ib1.le.nbc1 .and. lcbs(ib1).gt.0)
            is_core_boundary = is_core_boundary .or. ib.eq.lcbs(ib1)
            ib1 = ib1+1
          end do
          do is = 0, ns1-1
            if (bccon(is,ib).eq.1 .and. bcene(ib).eq.9) bcene(ib)=19
            if (bccon(is,ib).eq.1 .and. bceni(ib).eq.9) bceni(ib)=19
            if (bccon(is,ib).eq.21 .and.
     .          (streql(solps_version,'5.0') .or. 
     .           streql(solps_version,'5.1'))) then
              bccon(is,ib) = 24  ! must move number because of clash in numbering
            end if
          end do
          if (bceni(ib).eq.21 .and.
     .        (streql(solps_version,'5.0') .or. 
     .         streql(solps_version,'5.1'))) then
            bceni(ib) = 25  ! must move number because of clash in numbering
          end if
          if (bcpot(ib).eq.0) then
            if (bcchar(ib).eq.'N') then
              bcpot(ib) = 13 ! new treatment for automatic drift current contribution
            else if (bcchar(ib).eq.'S' .and. .not.is_core_boundary) then 
              bcpot(ib) = 13 ! new treatment for automatic drift current contribution
            else if (bcchar(ib).eq.'S' .and. is_core_boundary) then 
              bcpot(ib) = 12 ! new treatment for automatic drift current contribution
            end if
          end if  
        end do
        write(*,*) 'Writing ',trim(new_filename)
        write(nout(5),'(a)') trim(new_filename)
        open(35,file=new_filename)
        write(35,'(a)') '&BOUNDARY'
        write(35,'(1x,a,i2,a)') 'NBC=', nbc1,','
        line = 'BCCHAR= '
        len=8
        do ib = 1, nbc1
          line = line(1:len)//'  '''//bcchar(ib)//''','
          len=len+3+1+2
        end do
        write(35,'(1x,a)') line(1:len)
        write(hlp_frm,'(a,i2,a)') '(1x,a,',nbc1,'(1x,i4,'',''))'
        if(minval(bcpos1(1:nbc1)).gt.-2) then
          write(35, hlp_frm) 'BCPOS=  ',(bcpos1(ib)  ,ib=1,nbc1)
          write(35, hlp_frm) 'BCSTART=',(bcstart1(ib),ib=1,nbc1)
          write(35, hlp_frm) 'BCEND=  ',(bcend1(ib)  ,ib=1,nbc1)
        else
          do ib = 1, nbc1
            if(bcpos1(ib).gt.-2) then
              write(35,'(1x,a,i2,a,i4,a1)')
     .         'BCPOS(',ib,')=  ',bcpos1(ib),','
              write(35,'(1x,a,i2,a,i4,a1)')
     .         'BCSTART(',ib,')=',bcstart1(ib),','
              write(35,'(1x,a,i2,a,i4,a1)')
     .         'BCEND(',ib,')=  ',bcend1(ib),','
            else
              write(35,'(1x,a,i2,a,i4,a1)') 'BC_LIST_SIZE(',ib,')=',
     .                                       bc_list_size1(ib),','
              write(hlp_frm,'(a,i4,a)')
     .         '(1x,a,i2,a,',bc_list_size1(ib),'(1x,i4,'',''))'
              write(35, hlp_frm) 'BC_LIST_X(1,',ib,')=',
     .         bc_list_x1(1:bc_list_size1(ib),ib)
              write(35, hlp_frm) 'BC_LIST_Y(1,',ib,')=',
     .         bc_list_y1(1:bc_list_size1(ib),ib)
            end if
          end do
        endif
        write(hlp_frm,'(a,i2,a)') '(1x,a,',nbc1,'(1x,i4,'',''))'
        write(35, hlp_frm) 'BCENE=  ',(bcene(ib)  ,ib=1,nbc1)
        write(35, hlp_frm) 'BCENI=  ',(bceni(ib)  ,ib=1,nbc1)
        write(35, hlp_frm) 'BCPOT=  ',(bcpot(ib)  ,ib=1,nbc1)
        write(hlp_frm,'(a,i2,a)') '(1x,a,1p,',nbc1,'(1x,g10.3,'',''))'
        write(35, hlp_frm) 'ENEPAR(1,1)=',(enepar(ib,1),ib=1,nbc1)
        if (maxval(abs(enepar(1:nbc1,2))).gt.0.0_R8)
     >   write(35, hlp_frm) 'ENEPAR(1,2)=',(enepar(ib,2),ib=1,nbc1)
        write(35, hlp_frm) 'ENIPAR(1,1)=',(enipar(ib,1),ib=1,nbc1)
        if (maxval(abs(enipar(1:nbc1,2))).gt.0.0_R8)
     >   write(35, hlp_frm) 'ENIPAR(1,2)=',(enipar(ib,2),ib=1,nbc1)
        write(35, hlp_frm) 'POTPAR(1,1)=',(potpar(ib,1),ib=1,nbc1)
        if (maxval(abs(potpar(1:nbc1,2))).gt.0.0_R8)
     >   write(35, hlp_frm) 'POTPAR(1,2)=',(potpar(ib,2),ib=1,nbc1)
        write(hlp_frb,'(a,i3,a)') '(1x,a,i2,a,',ns1,'(1x,i3,'',''))'
        write(hlp_frp,'(a,i3,a)')
     .   '(1x,a,i2,a,1p,',ns1,'(1x,g10.3,'',''))'
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
          write(hlp_frm,'(a,i3,a)') '(1x,a,1p,',ns1,'(1x,g10.3,'',''))'
          write(35, hlp_frm) 'NIISO=',niiso(0:ns1-1)
          write(35, '(1x,a,1p,g10.3,a,g10.3,a,g10.3)')
     .     'TEISO=',teiso,', TIISO=',tiiso,', PHIISO=',phiiso
        endif
        if (boundary_time_mod.ne.0.0_R8) write(35,'(1x,a,f10.7,a)')
     >   'BOUNDARY_TIME_MOD=', boundary_time_mod,','
        if (boundary_time_switch.ne.0.0_R8) write(35,'(1x,a,f10.7,a)')
     >   'BOUNDARY_TIME_SWITCH=', boundary_time_switch,','
        if (boundary_time_mod.ne.0.0_R8) write(35,'(a,a,a)')
     >   'BOUNDARY_FILENAME=''', trim(boundary_filename),''','
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
        if (.not.streql(boundary_filename,'b2.boundary.parameters'))
     >   then
          new_filename =
     .     '../'//boundary_filename(1:lnblnk(boundary_filename))//'.2'
          filename = '../'//boundary_filename
          inquire(file=trim(filename),exist=todo)
          inquire(file=trim(new_filename),exist=found)
          if (found) write(6,*) trim(new_filename), ' already present!'
        else
          todo = .false.
          found = .true.
        endif
        boundary_namelist = 0
      enddo

*  ..write new b2.transport.parameters
      first_pass = .true.
      inquire(file='../b2.transport.parameters',exist=todo)
      new_filename = '../b2.transport.parameters.2'
      inquire(file=trim(new_filename),exist=found)
      if (found) write(6,*) trim(new_filename), ' already present!'
      do while((transport_namelist.eq.1.or.todo).and..not.found)
        if (transport_namelist.eq.0) call alloc_transport_namelist(ns)
        call read_b2mod_transport_namelist
        if (same_bundling) then
          if (first_pass) 
     >     write(nout(5),'(a)') 'b2.transport.parameters'
          if (.not.streql(transport_filename,'b2.transport.parameters'))
     >     write(nout(5),'(a)') trim(transport_filename)
        else
          write(nout(5),'(a)') trim(new_filename)
          write(*,*) 'Writing ',trim(new_filename)
          write(hlp_frm,'(a,i3,a)')
     .     '(1x,a,i2,a,a,1p,',ns1,'(g10.3,'','',1x))'
          open(99,file=trim(new_filename))
          write(99,'(a10)') '&TRANSPORT'
          if(flag_dna.ne.0) then
            cftmp(0:ns-1)=parm_dna(0:ns-1)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            parm_dna(0:ns1-1)=cftmp(0:ns1-1)
            write(99,hlp_frm)
     .       'FLAG_DNA=',flag_dna,',',' PARM_DNA=',parm_dna(0:ns1-1)
          endif
          if(flag_dpa.ne.0) then
            cftmp(0:ns-1)=parm_dpa(0:ns-1)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          parm_dpa(0:ns1-1)=cftmp(0:ns1-1)
          write(99,hlp_frm)
     .     'FLAG_DPA=',flag_dpa,',',' PARM_DPA=',parm_dpa(0:ns1-1)
          endif
          if(flag_vla.ne.0) then
            cftmp(0:ns-1)=parm_vla(0:ns-1)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            parm_vla(0:ns1-1)=cftmp(0:ns1-1)
            write(99,hlp_frm)
     .       'FLAG_VLA=',flag_vla,',',' PARM_VLA=',parm_vla(0:ns1-1)
          endif
          if(flag_vsa.ne.0) then
            cftmp(0:ns-1)=parm_vsa(0:ns-1)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            parm_vsa(0:ns1-1)=cftmp(0:ns1-1)
            write(99,hlp_frm)
     .       'FLAG_VSA=',flag_vsa,',',' PARM_VSA=',parm_vsa(0:ns1-1)
          endif
          if(flag_hci.ne.0) then
            cftmp(0:ns-1)=parm_hci(0:ns-1)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            parm_hci(0:ns1-1)=cftmp(0:ns1-1)
            write(99,hlp_frm)
     .       'FLAG_HCI=',flag_hci,',',' PARM_HCI=',parm_hci(0:ns1-1)
          endif
          if(flag_hce.ne.0) write(99,'(1x,a,i2,a,a,1p,1g10.3,a)')
     .     'FLAG_HCE=',flag_hce,',',' PARM_HCE=',parm_hce,','
          if(flag_sig.ne.0) write(99,'(1x,a,i2,a,a,1p,1g10.3,a)')
     .     'FLAG_SIG=',flag_sig,',',' PARM_SIG=',parm_sig,','
          if(flag_alf.ne.0) write(99,'(1x,a,i2,a,a,1p,1g10.3,a)')
     .     'FLAG_ALF=',flag_alf,',',' PARM_ALF=',parm_alf,','
          if(cflme.ge.0.0_R8) write(99,'(1x,a,1p,1g10.3,a)')
     .     'CFLME=',cflme,','
          if(cflmi.ge.0.0_R8) write(99,'(1x,a,1p,1g10.3,a)')
     .     'CFLMI=',cflmi,','
          if(cflmv.ge.0.0_R8) write(99,'(1x,a,1p,1g10.3,a)')
     .     'CFLMV=',cflmv,','
          if(transport_time_mod.gt.0.0_R8) then
            write(99,'(1x,a,a,a)')
     .       'TRANSPORT_FILENAME=''',trim(transport_filename),''','
            write(99,'(1x,a,1p,1g10.3,a)')
     .       'TRANSPORT_TIME_MOD=',transport_time_mod,','
            write(99,'(1x,a,1p,1g10.3,a)')
     .       'TRANSPORT_TIME_SWITCH=',transport_time_switch,','
          endif
          if(abs(maxval(vout_cnv(0:ns-1))).gt.0.0_R8) then
            write(hlp_frm,'(a,i3,a)')
     .       '(1x,a,1p,',ns1,'(1x,g10.3,'',''))'
            cftmp(0:ns-1)=vout_cnv(0:ns-1)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            vout_cnv(0:ns1-1)=cftmp(0:ns1-1)
            write(99,hlp_frm)
     .       'VOUT_CNV=',vout_cnv(0:ns1-1)
            cftmp(0:ns-1)=pw0_cnv(0:ns-1)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            pw0_cnv(0:ns1-1)=cftmp(0:ns1-1)
            write(99,hlp_frm)
     .       'PW0_CNV=',pw0_cnv(0:ns1-1)
            cftmp(0:ns-1)=pw1_cnv(0:ns-1)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            pw1_cnv(0:ns1-1)=cftmp(0:ns1-1)
            write(99,hlp_frm)
     .       'PW1_CNV=',pw1_cnv(0:ns1-1)
            cftmp(0:ns-1)=pw2_cnv(0:ns-1)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            pw2_cnv(0:ns1-1)=cftmp(0:ns1-1)
            write(99,hlp_frm)
     .       'PW2_CNV=',pw2_cnv(0:ns1-1)
          endif
          write(99,'(a1)') '/'
          close(99)
        end if
        if (.not.streql(transport_filename,'b2.transport.parameters'))
     >  then
          first_pass = .false.
          new_filename =
     .     '../'//transport_filename(1:lnblnk(transport_filename))//'.2'
          filename='../'//transport_filename
          inquire(file=trim(filename),exist=todo)
          inquire(file=trim(new_filename),exist=found)
          if (found) write(6,*) trim(new_filename), ' already present!'
        else
          todo = .false.
          found = .true.
        endif
        transport_namelist = 0
      enddo

*  ..write new b2.transport.inputfile
      first_pass = .true.
      filename='../b2.transport.inputfile'
      inquire(file=trim(filename),exist=todo)
      new_filename='../b2.transport.inputfile.2'
      inquire(file=trim(new_filename),exist=found)
      if (found) write(6,*) trim(new_filename), ' already present!'
      do while((transport_inputfile.eq.1.or.todo).and..not.found)
        if (.not.allocated(px)) then
          allocate(px(1:ny1+2),pf(1:ny1+2))
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
          transport_filename='b2.transport.inputfile'
          transport_time_mod=0.0_R8
          transport_time_switch=0.0_R8
          allocate(n1data(1:nkind_data,1:nkind_coeff,0:nsm-1))
          allocate(t1data(3,1:nym+2,1:nkind_coeff,0:nsm-1))
          allocate(ad1spec(ns1,1:nkind_coeff,0:ns1-1))
          allocate(pelm(1:ny1+2))
          n1data(:,:,:)=0
          t1data(:,:,:,:)=0.0_R8
          ad1spec(:,:,:)=-5
        end if
        open(99,file=trim(filename))
        read(99,transport)
        close(99)
        if(same_bundling.and.nx.eq.nx1.and.ny.eq.ny1) then
          if (first_pass) write(nout(5),'(a)') 'b2.transport.inputfile'
          if (.not.streql(transport_filename,'b2.transport.inputfile'))
     >     write(nout(5),'(a)') trim(transport_filename)
        else
          write(nout(5),'(a)') trim(new_filename)
          write(*,*) 'Writing ',trim(new_filename)
          open(99,file=trim(new_filename))
          write(99,'(a10)') '&TRANSPORT'
          do j = 1, nkind_coeff
            do i = 1, nkind_data
cxpb  We do the mesh size transformation first
              do is = 0, ns-1
                prdone = .false.
                ndat = ndata(i,j,is)
                if (ndat.gt.0 .and. is1(is).ge.0 .and. is1(is).lt.ns1)
     .           then
                  if (ndat.gt.ny1+2) then  ! need to coarsen the profile according to the new grid
                    if (.not.prdone) then
                      allocate(prdsa(1:ny1+2),pra(1:ny1+2))
                      allocate(r(1:ny+2),f(1:ny+2),felm(1:ny+2))
                      allocate(dprof(1:ny+2),delm(1:ny+2))
                      prdsa(1)=sqrt((cr1(jxa1,-1)-
     1                 0.5_R8*(crx1(jxa1,-1,0)+crx1(jxa1,-1,1)))**2+
     2                  (cz1(jxa1,-1)-
     3                 0.5_R8*(cry1(jxa1,-1,0)+cry1(jxa1,-1,1)))**2)
                      do iy = 0, ny1
                        prdsa(iy+2)=prdsa(iy+1)+
     1                   sqrt((cr1(jxa1,iy)-cr1(jxa1,iy-1))**2+
     2                        (cz1(jxa1,iy)-cz1(jxa1,iy-1))**2)
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
                    call E01BEF(ndat,r,f,dprof,ifail)
                    call E01BEF(ndat,r,felm,delm,ifail)
                    nfitanf=1
                    do iy=1,ny1+2
                      pra(iy)=prdsa(iy)
                      if(pra(iy).lt.r(1)) nfitanf=nfitanf+1
                      nfitend=iy-1
                      if (pra(iy).GT.r(ndat)) exit
                      nfitend=iy
                    enddo
                    call xertst(nfitanf.lt.nfitend,
     .               'Could not complete transport profile!')
                    ndim=nfitend+1-nfitanf
                    do iy = 1, ndim
                      px(i) = pra(nfitanf-1+iy)
                    end do
                    call E01BFF(ndat,r,f,dprof,ndim,px,pf,ifail)
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
                  write(99,'(1x,a,i1,a,i1,a,i3,a,i3,a)')
     .             'NDATA(',i,',',j,',',is1(is),')=',
     .              n1data(i,j,is1(is)),','
                  iy = 0
cxpb Check if we need to include new additional species because of de-bundling of is
                  if (bundle_frac(is1(is),is).lt.1.0_R8) then
                    do js = 0, ns1-1
                      if (js.ne.is1(is) .and.
     &                    bundle_frac(js,is).ne.0.0_R8) then
                        iy = iy + 1
                        ad1spec(iy,j,is1(is)) = js
                        write(99,'(1x,a,i3,a,i1,a,i3,a,i3,a)')
     .                   'ADDSPEC(',iy,',',j,',',is1(is),')=',
     .                    ad1spec(iy,j,is1(is)),','
                      endif
                    enddo
                  endif
                  do ix = 1, ns
cxpb Convert old additional species
                    if (addspec(ix,j,is).ge.0) then
                      todo = .true.
                      ib = 1
                      do while (ib.le.ns1)
                        if (ad1spec(ib,j,is1(is)).ne.-5) then
                          ib = ib + 1
                        else
                          iy = 1
                          do while (iy.lt.ns1)
cxpb We must ensure we have not already listed this additional species
                            found = .false.
                            iss = 1
                            do while (.not.found .and. iss.lt.ib)
                              found = ad1spec(iss,j,is1(is)).eq.iy
                              iss = iss + 1
                            end do
                            if(
     .                       bundle_frac(iy,addspec(ix,j,is)).gt.0.0_R8
     .                       .and. .not.found) then
                              ad1spec(ib,j,is1(is)) = iy
                              write(99,'(1x,a,i3,a,i1,a,i3,a,i3,a)')
     .                         'ADDSPEC(',ib,',',j,',',is1(is),')=',
     .                          ad1spec(ib,j,is1(is)),','
                              ib = ib + 1
                            endif
                            iy = iy + 1
                          enddo
                          ib = ns1 + 1
                        endif
                      enddo
                    endif
                  enddo
                  do ix = 1, ndim
                    write(99,'(1x,a,i3,a,i1,a,i3,a,1p,3(g10.3,a))')
     .               'TDATA(1,',ix,',',j,',',is1(is),')=',
     .               t1data(1,ix,j,is1(is)),',',
     .               t1data(2,ix,j,is1(is)),',',
     .               t1data(3,ix,j,is1(is)),','
                  enddo
                endif
                if (prdone) deallocate(felm,delm,prdsa,r,f,pra,dprof)
              enddo
            enddo
            todo = .false.
            do ireg = 1, nnreg(0)
              todo = todo .or. .not.region_flags(ireg,j)
            enddo
            write(hlp_frm, '(a,i1,a)') '(1x,a,i1,a,',nnreg(0),'l1,a)'
            if (todo) write(99,hlp_frm)
     .       'REGION_FLAGS(1,',j,')=',
     .       (region_flags(ireg,j),ireg=1,nnreg(0)),','
          enddo
          write(99,'(1x,a,l1,a)') 'NO_PFLUX=',no_pflux,','
          if(transport_time_mod.gt.0.0_R8) then
            write(99,'(1x,a,a,a)')
     .       'TRANSPORT_FILENAME=''',trim(transport_filename),''','
            write(99,'(1x,a,1p,1g10.3,a)')
     .       'TRANSPORT_TIME_MOD=',transport_time_mod,','
            write(99,'(1x,a,1p,1g10.3,a)')
     .       'TRANSPORT_TIME_SWITCH=',transport_time_switch,','
          endif
          if (poloidal_scaling) write(99,'(1x,a,l1,a)')
     .     'POLOIDAL_SCALING=',poloidal_scaling,','
          if (scaling_ix_begin.ne.-2) then
            scaling_ix_begin = conformal_convert(non_commensurate,1,
     .       scaling_ix_begin,jsep,nx,ny,nx1,ny1,1)
            write(99,'(1x,a,i4,a)')
     .       'SCALING_IX_BEGIN=',scaling_ix_begin,','
          endif
          if (scaling_ix_end.ne.-2) then
            scaling_ix_end = conformal_convert(non_commensurate,1,
     .       scaling_ix_end,jsep,nx,ny,nx1,ny1,3)
            write(99,'(1x,a,i4,a)')
     .       'SCALING_IX_END=',scaling_ix_end,','
          endif
          if (scaling_strength.ne.0.0_R8) write(99,'(1x,a,1p,g10.3,a)')
     .     'SCALING_STRENGTH=',scaling_strength,','
          if (scaling_width.ne.0.0_R8) write(99,'(1x,a,1p,g10.3,a)')
     .     'SCALING_WIDTH=',scaling_width,','
          if (elm_time_begin.ne.0.0_R8) write(99,'(1x,a,1p,g10.3,a)')
     .     'ELM_TIME_BEGIN=',elm_time_begin,','
          if (elm_time_end.ne.0.0_R8) write(99,'(1x,a,1p,g10.3,a)')
     .     'ELM_TIME_END=',elm_time_end,','
          if (elm_time_period.ne.0.0_R8) write(99,'(1x,a,1p,g10.3,a)')
     .     'ELM_TIME_PERIOD=',elm_time_period,','
          if (elm_ix_begin.ne.-2) then
            elm_ix_begin = conformal_convert(non_commensurate,1,
     .       elm_ix_begin,jsep,nx,ny,nx1,ny1,1)
            write(99,'(1x,a,i4,a)')
     .       'ELM_IX_BEGIN=',elm_ix_begin,','
          endif
          if (elm_ix_end.ne.-2) then
            elm_ix_end = conformal_convert(non_commensurate,1,
     .       elm_ix_end,jsep,nx,ny,nx1,ny1,3)
            write(99,'(1x,a,i4,a)')
     .       'ELM_IX_END=',elm_ix_end,','
          endif
          write(99,'(a1)') '/'
          close(99)
        end if
        if (.not.streql(transport_filename,'b2.transport.inputfile'))
     >   then
          first_pass = .false.
          new_filename =
     .     '../'//transport_filename(1:lnblnk(transport_filename))//'.2'
          filename='../'//transport_filename
          inquire(file=trim(filename),exist=todo)
          inquire(file=trim(new_filename),exist=found)
          if (found) write(6,*) trim(new_filename), ' already present!'
        else
          todo = .false.
          found = .true.
        endif
        transport_inputfile = 0
        deallocate(t1data)
        deallocate(n1data)
        deallocate(ad1spec)
        deallocate(pelm)
        deallocate(px)
        deallocate(pf)
      enddo

*  ..write new b2.transport_models_save.parameters
      inquire(file='../b2.transport_models_save.parameters',exist=todo)
      if (todo .and. (ny.ne.ny1 .or. non_commensurate.eq.1)) then
        call alloc_b2mod_transport_models(nx,ny,nsm)
        call read_b2mod_transport_models(99,'    ',todo)
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
          write(nout(5),'(a)') 'b2.transport_models_save.parameters.2'
          open(99,file='../b2.transport_models_save.parameters.2')
          write(99,'(a9)') '&TRANSPORT_MODELS'
          write(hlp_frm,'(a,i4,a)')
     .     '(1x,a,1p,',ny1+2,'(1x,g10.3,'',''))'
          write(99,hlp_frm)
     .     'ETA_HCE_MULT=',eta_hce_mult1(-1:ny1)
          write(99,'(a1)') '/'
          close(99)
        endif
        deallocate(dsa,dsa1,deta)
        deallocate(eta_hce_mult1)
      endif

*  ..write new b2.sources.profile
      axdone = .false.
      nsdata(:,:,:)=0
      sdata(:,:,:,:)=0.0_R8
      nxdata(:,:,:)=0
      xdata(:,:,:,:)=0.0_R8
      divheat=0.0_R8
      sources_filename='b2.sources.profile'
      sources_time_mod=0.0_R8
      sources_time_switch=0.0_R8
      allocate(ns1data(1:nkind_data,1:nkind_source,0:nsm-1))
      allocate(s1data(2,1:DEF_NYD+2,1:nkind_source,0:nsm-1))
      allocate(nx1data(1:nkind_data,1:nkind_source,0:nsm-1))
      allocate(x1data(2,1:DEF_NXD+2,1:nkind_source,0:nsm-1))
      allocate(px(1:ny1+2),pf(1:ny1+2))
      ns1data(:,:,:)=0
      s1data(:,:,:,:)=0.0_R8
      nx1data(:,:,:)=0
      x1data(:,:,:,:)=0.0_R8
      inquire(file='../b2.sources.profile',exist=todo)
      new_filename = '../b2.sources.profile.2'
      inquire(file=trim(new_filename),exist=found)
      if (found) write(6,*) trim(new_filename), ' already present!'
      do while((sources_inputfile.eq.1.or.todo).and..not.found)
        filename = '../'//sources_filename
        open(92,file=trim(filename))
        read(92,profile)
        close(92)
        write(*,*) 'Writing ',trim(new_filename)
        write(nout(5),'(a)') trim(new_filename)
        open(99,file=new_filename)
        write(99,'(a9)') '&PROFILE'
        do ikind = 1, nkind_data
          do j = 1, nkind_source
            do is = 0, ns-1
cxpb First we do the radial profile
              prdone=.false.
              ndat=nsdata(ikind,j,is)
              if (ndat.gt.0.and.is1(is).ge.0.and.is1(is).lt.ns1) then
                clash(1:ns-1) = .false.
                shared(0:ns-1)= .false.
                w1d(0)=bundle_frac(is1(is),is)
                ks = 0
                do js = 0, ns-1
                  if (js.ne.is .and. nsdata(ikind,j,js).gt.0 .and.
     &                bundle_frac(is1(is),js).ne.0.0_R8) then  ! must check the profiles do not clash
                    shared(js)=.true.
                    ks = ks + 1
                    w1d(ks)=bundle_frac(is1(is),js)
cxpb Do we have the same number of points?
                    nds(ks) = nsdata(ikind,j,js)
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
                  allocate(spr(1:nnc*(ny+2+2)))
                  ndim = ndat
                  allocate(s1d(2,1:ndim,1:nkind_source,0:nsm-1))
                  spr(1:ndat) = sr1(1:ndat)
                  do js = 0, ns-1
                    if (shared(js)) then
                      ik = ik+1
                      if (clash(ik)) then
                        allocate(sr2(1:nds(ik)),sf2(1:nds(ik)))
                        allocate(sd2(1:nds(ik)))
                        allocate(spf2(1:ndim))
                        sr2(1:nds(ik)) = s1data(1,1:nds(ik),j,js)
                        sf2(1:nds(ik)) = s1data(2,1:nds(ik),j,js)
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
                  call xertst(nfitanf.lt.nfitend,
     .             'Could not complete initial source profile!')
                  ndm2=nfitend+1-nfitanf
                  call E01BFF(ndat,sr1,sf1,sd1,ndm2,
     .             spr(nfitanf:nfitend),spf1(nfitanf:nfitend),ifail)
                  spf1(1:nfitanf-1)=0.0_R8
                  spf1(nfitend+1:ndim)=0.0_R8
                  s1d(1,1:ndim,j,is)=spr(1:ndim)
                  s1d(2,1:ndim,j,is)=spf1(1:ndim)*w1d(0)
                  deallocate(sr1,sf1,sd1,spf1)
                  ns1data(ikind,j,is)=ndim
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
                        call xertst(nfitanf.lt.nfitend,
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
                  ns1data(ikind,j,is)=ndat
                  s1d(1,1:ndat,j,is)=s1data(1,1:ndat,j,is)
                  s1d(2,1:ndat,j,is)=s1data(2,1:ndat,j,is)
                endif
                do js = 0, ns-1
                  if (bundle_frac(is1(is),js).eq.1.0_R8 .and. is.ne.js)
     >             then
                    ns1data(ikind,j,js)=0   ! This profile has been absorbed in the current one
                  else if (bundle_frac(is1(is),js).gt.0.0_R8) then
                    if (ndat.gt.ny1+2) then  ! need to coarsen the profile according to the new grid
                      if (.not.prdone) then
                        allocate(prdsa(1:ny1+2),pra(1:ny1+2))
                        allocate(r(1:ny+2),f(1:ny+2))
                        allocate(ds(1:ny+2))
                        prdsa(1)=sqrt((cr1(ixref1,-1)-
     1                   0.5_R8*(crx1(ixref1,-1,0)+
     .                           crx1(ixref1,-1,1)))**2+
     2                    (cz1(ixref1,-1)-
     3                   0.5_R8*(cry1(ixref1,-1,0)+
     .                           cry1(ixref1,-1,1)))**2)
                        do iy = 0, ny1
                          prdsa(iy+2)=prdsa(iy+1)+
     1                     sqrt((cr1(ixref1,iy)-cr1(ixref1,iy-1))**2+
     2                          (cz1(ixref1,iy)-cz1(ixref1,iy-1))**2)
                        enddo
                        ds_offset=
     .                   (prdsa(iyref1+2)+prdsa(iyref1+3))/2.0_R8
                        do iy = 1, ny1+2
                          prdsa(iy) = prdsa(iy) - ds_offset
                        enddo
                        prdone = .true.
                      endif
                      r(1:ndat) = s1d(1,1:ndat,j,is)
                      f(1:ndat) = s1d(2,1:ndat,j,is)
                      call E01BEF(ndat,r,f,ds,ifail)
                      nfitanf=1
                      do iy=1,ny1+2
                        pra(iy)=prdsa(iy)
                        if(pra(iy).lt.r(1)) nfitanf=nfitanf+1
                        nfitend=iy-1
                        if (pra(iy).GT.r(ndat)) exit
                        nfitend=iy
                      enddo
                      call xertst(nfitanf.lt.nfitend,
     .                 'Could not complete rescaled source profile!')
                      ndim=nfitend+1-nfitanf
                      do iy = 1, ndim
                        px(iy) = pra(nfitanf-1+iy)
                      end do
                      call E01BFF(ndat,r,f,ds,ndim,px,pf,ifail)
                    else  ! simple case
                      ndim = ndat
                      px(1:ndim) = s1d(1,1:ndim,j,is)
                      pf(1:ndim) = s1d(2,1:ndim,j,is)
                    endif
                    ns1data(ikind,j,is1(is))=ndim
                    s1data(1,1:ndim,j,is1(is))=px(1:ndim)  ! profile position
                    s1data(2,1:ndim,j,is1(is))=pf(1:ndim)  ! profile value
                  endif
                enddo
                deallocate(s1d)
                if (prdone) deallocate(prdsa,pra,r,f,ds)
              endif
              if (ns1data(ikind,j,is1(is)).gt.0) then
                write(99,'(1x,a,i1,a,i1,a,i3,a,i3,a)')
     .           'NSDATA(',ikind,',',j,',',is1(is),')=',
     .            ns1data(ikind,j,is1(is)),','
                do ix = 1, ns1data(ikind,j,is1(is))
                  write(99,'(1x,a,i3,a,i1,a,i3,a,1p,2(g12.5,a))')
     .             'SDATA(1,',ix,',',j,',',is1(is),')=',
     .             s1data(1,ix,j,is1(is)),',',
     .             s1data(2,ix,j,is1(is)),','
                enddo
                write(99,1028)
              endif
cxpb Second we do the axial profile
              ndat=nxdata(ikind,j,is)
              axdone=.false.
              if (ndat.gt.0.and.is1(is).ge.0.and.is1(is).lt.ns1) then
                clash(1:ns-1) = .false.
                shared(0:ns-1)= .false.
                w1d(0)=bundle_frac(is1(is),is)
                ks = 0
                do js = 0, ns-1
                  if (js.ne.is .and. nxdata(ikind,j,js).gt.0 .and.
     &                bundle_frac(is1(is),js).ne.0.0_R8) then  ! must check the profiles do not clash
                    shared(js)=.true.
                    ks = ks + 1
                    w1d(ks)=bundle_frac(is1(is),js)
cxpb Do we have the same number of points?
                    nds(ks) = nxdata(ikind,j,js)
                    clash(ks) = clash(ks) .or. ndat.ne.nds(ks)
                    if (.not.clash(ks)) then
cxpb Are the points in the same place?
                      do ix = 1, ndat
                        clash(ks) = clash(ks) .or.
     &                   abs(xdata(1,ix,j,is)-xdata(1,ix,j,js)).le.
     .                   1.0e-6_R8
                      enddo
                    endif
                  endif
                enddo
                ik = 0
                inc = 0
                x1data(1,1:ndat,j,is)=xdata(1,1:ndat,j,is)
                x1data(2,1:ndat,j,is)=xdata(2,1:ndat,j,is)*w1d(0)
                do js = 0, ns-1
                  if (shared(js)) then

                    ik = ik+1
                    if (.not.clash(ik)) then
cxpb We can sum up all the profiles that do not clash
                      inc = inc+1
                      x1data(2,1:ndat,j,is)=xdata(2,1:ndat,j,is)+
     +                                      xdata(2,1:ndat,j,js)*w1d(ik)
                    endif
                  endif
                enddo
cxpb Some profiles clash: we must compute a new combined profile
                if (inc.lt.ks) then
                  nnc = ks-inc
                  ik = 0
                  allocate(xr1(1:(nx+2)),xf1(1:(nx+2)))
                  allocate(xd1(1:(nx+2)))
                  xr1(1:ndat) = xdata(1,1:ndat,j,is)
                  xf1(1:ndat) = xdata(2,1:ndat,j,is)
                  call E01BEF(ndat,xr1,xf1,xd1,ifail)
cxpb Calculate axial locations of combined axial profile xpr(1:ndim)
cxpb This will need to be checked !!!
                  allocate(xpr(1:nnc*(nx+2)))
                  ndim = ndat
                  xpr(1:ndat) = xr1(1:ndat)
                  do js = 0, ns-1
                    if (shared(js)) then
                      ik = ik+1
                      if (clash(ik)) then
                        allocate(xr2(1:nds(ik)),xf2(1:nds(ik)))
                        allocate(xd2(1:nds(ik)))
                        allocate(xpf2(1:ndim))
                        xr2(1:nds(ik)) = x1data(1,1:nds(ik),j,js)
                        xf2(1:nds(ik)) = x1data(2,1:nds(ik),j,js)
                        i1 = 1
                        i2 = 1
                        do while(i2.le.nds(ik))
                          if (abs(xpr(i1)-xr2(i2)).lt.1.0e-6_R8) then  ! shared point
                            i1 = i1+1
                            i2 = i2+1
                          else if (xpr(i1).lt.xr2(i2)) then  ! catch up with 2nd profile
                            i1 = i1+1
                          else ! add a new point to the profile
                            do ix=ndim,i1,-1
                              xpr(ix+1)=xpr(ix)
                            enddo
                            xpr(i1)=xr2(i2)
                            ndim=ndim+1
                            i2 = i2+1
                          endif
                        enddo
                      endif
                    endif
                  enddo
                  allocate(xpf1(1:ndim))
                  nfitanf=1
                  do ix = 1, ndim
                    if(xpr(ix).lt.xr1(1)) nfitanf=nfitanf+1
                    nfitend=iy-1
                    if(xpr(ix).gt.xr1(ndat)) exit
                    nfitend=ix
                  enddo
                  call xertst(nfitanf.lt.nfitend,
     .             'Could not complete initial axial profile!')
                  ndm2=nfitend+1-nfitanf
                  call E01BFF(ndat,xr1,xf1,xd1,ndm2,
     .             xpr(nfitanf:nfitend),xpf1(nfitanf:nfitend),ifail)
                  xpf1(1:nfitanf-1)=0.0_R8
                  xpf1(nfitend+1:ndim)=0.0_R8
                  allocate(x1d(2,1:ndim,1:nkind_source,0:nsm-1))
                  x1d(1,1:ndim,j,is)=xpr(1:ndim)
                  x1d(2,1:ndim,j,is)=xpf1(1:ndim)*w1d(0)
                  deallocate(xr1,xf1,xd1,xpf1)
                  nx1data(ikind,j,is)=ndim
                  ndat=ndim
                  ik = 0
                  do js = 0, ns-1
                    if (shared(js)) then
                      ik = ik+1
                      if (clash(ik)) then
                        allocate(xr2(1:nds(ik)),xf2(1:nds(ik)))
                        allocate(xd2(1:nds(ik)))
                        allocate(xpf2(1:ndim))
                        xr2(1:nds(ik)) = x1data(1,1:nds(ik),j,js)
                        xf2(1:nds(ik)) = x1data(2,1:nds(ik),j,js)
                        nfitanf=1
                        do ix = 1, ndim
                          if(xpr(ix).lt.xr2(1)) nfitanf=nfitanf+1
                          nfitend=ix-1
                          if(xpr(ix).gt.xr2(nds(ik))) exit
                          nfitend=ix
                        enddo
                        call xertst(nfitanf.lt.nfitend,
     .                   'Could not complete additional axial profile!')
                        ndm2=nfitend+1-nfitanf
                        call E01BEF(nds(ik),xr2,xf2,xd2,ifail)
                        call E01BFF(nds(ik),xr2,xf2,xd2,ndm2,
     .                   xpr(nfitanf:nfitend),xpf2(nfitanf:nfitend),
     .                   ifail)
                        xpf2(1:nfitanf-1)=0.0_R8
                        xpf2(nfitend+1:ndim)=0.0_R8
                        x1d(2,1:ndim,j,is)=
     .                   x1d(2,1:ndim,j,is)+xpf2(1:ndim)*w1d(ik)
                        deallocate(xr2,xf2,xd2,xpf2)
                      endif
                    endif
                  enddo
                  deallocate(xpr)
                else  ! no clash : copy over data
                  allocate(x1d(2,1:ndat,1:nkind_source,0:nsm-1))
                  nx1data(ikind,j,is)=ndat
                  x1d(1,1:ndat,j,is)=x1data(1,1:ndat,j,is)
                  x1d(2,1:ndat,j,is)=x1data(2,1:ndat,j,is)
                endif
                do js = 0, ns-1
                  if (bundle_frac(is1(is),js).eq.1.0_R8 .and. is.ne.js)
     >             then
                    nx1data(ikind,j,js)=0   ! This profile has been absorbed in the current one
                  else if (bundle_frac(is1(is),js).gt.0.0_R8) then
                    if (ndat.gt.nx1+2) then  ! need to coarsen the profile according to the new grid
                      if (.not.axdone) then
                        allocate(pxa(-1:nx1))
                        allocate(l_conn(-1:nx1))
                        allocate(x(1:nx+2),fx(1:nx+2))
                        allocate(dx(1:nx+2))
                        ireg = mod(region1(ixref1,iyref1,0),4)
                        if (ikind.eq.1) then
                          l_total = 0.0_R8
                          pxa(-1) = 0.0_R8
                          l_conn(-1) = 0.0_R8
                          do i = 0, nx1-1
                            l_conn(i) = hx1(i,iyref1)/
     1                       (sqrt(bb1(i,iyref1,0)**2+
     2                             bb1(i,iyref1,1)**2)/bb1(i,iyref1,3))
                            pxa(i) =
     1                       pxa(i-1) + 0.5_R8*(l_conn(i)+l_conn(i-1))
                            l_total = l_total + l_conn(i)
                          end do
                          pxa(nx1) = pxa(nx1-1) + 0.5_R8*l_conn(nx1-1)
                          l_ref = pxa(ixref1)
                          if (ireg.eq.1 .and. nnreg(0).ge.4) then
                            do i = -1, nx1
                              if (mod(region1(i,iyref1,0),4).eq.ireg)
     >                         then
                                pxa(i) = (pxa(i) - l_ref) / l_total
                                if (pxa(i).lt.0.0_R8)
     >                           pxa(i) = pxa(i) + 1.0_R8
                              else
                                pxa(i) = 0.0_R8
                              end if
                            end do
                          else
                            do i = -1, nx1
                              pxa(i) = pxa(i) / l_total
                            end do
                          end if
                        else if (ikind.eq.2) then
! Establish profile in computational space
                          pxa(-1) = 0.0_R8
                          do i = 0, nx1-1
                            pxa(i) = (i + 0.5_R8) / (real(nx1))
                          end do
                          pxa(nx1) = 1.0_R8
                        end if
                        nfitanf = -1
                        do i = 0, nx1
                          if (pxa(i).lt.x(1)) nfitanf = nfitanf+1
                          nfitend=i-1
                          if (pxa(i).gt.x(ndat)) exit
                          nfitend=i
                        end do
                        ndim=nfitend+1-nfitanf
                        axdone = .true.
                      endif
                      xx(1:ndat) = x1d(1,1:ndat,j,is)
                      fx(1:ndat) = x1d(2,1:ndat,j,is)
                      call E01BEF(ndat,xx,fx,dx,ifail)
                      nfitanf=1
                      do ix=1,nx1+2
                        pxa(ix)=pxdsa(ix)
                        if(pxa(ix).lt.x(1)) nfitanf=nfitanf+1
                        nfitend=ix-1
                        if (pxa(ix).GT.x(ndat)) exit
                        nfitend=ix
                      enddo
                      call xertst(nfitanf.lt.nfitend,
     .                 'Could not complete rescaled axial profile!')
                      ndim=nfitend+1-nfitanf
                      do ix = 1, ndim
                        px(ix) = pxa(nfitanf-1+ix)
                      end do
                      call E01BFF(ndat,xx,fx,dx,ndim,px,pf,ifail)
                      nx1data(ikind,j,is1(is))=ndim
                      x1data(1,1:ndim,j,is1(is)) = px(1:ndim)  ! profile position
                      x1data(2,1:ndim,j,is1(is)) = pf(1:ndim)  ! profile value
                    else  ! simple case
                      nx1data(ikind,j,is1(is))=ndat
                      x1data(1,1:ndat,j,is1(is)) = x1d(1,1:ndat,j,is)  ! profile position
                      x1data(2,1:ndat,j,is1(is)) = x1d(2,1:ndat,j,is)  ! profile value
                    endif
                    if (axdone) deallocate(pxa,x,fx,dx)
                  endif
                enddo
                deallocate(x1d)
              endif
              if (nx1data(ikind,j,is1(is)).gt.0) then
                write(99,'(1x,a,i1,a,i1,a,i3,a,i3,a)')
     .           'NXDATA(',ikind,',',j,',',is1(is),')=',
     .            nx1data(ikind,j,is1(is)),','
                do ix = 1, nx1data(ikind,j,is1(is))
                  write(99,'(1x,a,i3,a,i1,a,i3,a,1p,2(g12.5,a))')
     .             'XDATA(1,',ix,',',j,',',is1(is),')=',
     .              x1data(1,ix,j,is1(is)),',',
     .              x1data(2,ix,j,is1(is)),','
                enddo
                write(99,1028)
              endif
            enddo
          enddo
        enddo
        if (divheat.ne.0.0_R8) write(99,'(1x,a,1p,1g10.3,a)')
     .   'DIVHEAT=',divheat,','
        if(sources_time_mod.gt.0.0_R8) then
          write(99,'(1x,a,1p,1g10.3,a)')
     .     'SOURCES_TIME_MOD=',sources_time_mod,','
          write(99,'(1x,a,1p,1g10.3,a)')
     .     'SOURCES_TIME_SWITCH=',sources_time_switch,','
          write(99,'(1x,a,a,a)')
     .     'SOURCES_FILENAME=''',trim(sources_filename),''','
        endif
        write(99,'(a1)') '/'
        close(99)
        if (.not.streql(sources_filename,'b2.sources.profile'))
     >   then
          new_filename =
     .     '../'//sources_filename(1:lnblnk(sources_filename))//'.2'
          filename = '../'//sources_filename
          inquire(file=trim(filename),exist=todo)
          inquire(file=trim(new_filename),exist=found)
          if (found) write(6,*) trim(new_filename), ' already present!'
        else
          todo = .false.
          found = .true.
        endif
        sources_inputfile = 0
      enddo
      deallocate(s1data)
      deallocate(ns1data)
      deallocate(x1data)
      deallocate(nx1data)
      deallocate(px,pf)

*  ..write new b2.numerics.parameters
      first_pass = .true.
      inquire(file='../b2.numerics.parameters',exist=todo)
      new_filename = '../b2.numerics.parameters.2'
      inquire(file=new_filename,exist=found)
      if (found) write(6,*) trim(new_filename), ' already present!'
      do while((numerics_namelist.ge.1.or.todo).and..not.found)
        if (numerics_namelist.eq.0) call alloc_b2mod_numerics(nx,ny,ns)
        call read_b2mod_numerics_namelist(nx,ny)
        if (same_bundling) then
          if (first_pass) 
     >     write(nout(5),'(a)') 'b2.numerics.parameters'
          if (.not.streql(numerics_filename,'b2.numerics.parameters'))
     >     write(nout(5),'(a)') trim(numerics_filename)
        else
          write(*,*) 'Writing ',trim(new_filename)
          open(99,file=new_filename)
          write(nout(5),'(a)') trim(new_filename)
          write(99,'(a9)') '&NUMERICS'
          write(hlp_frm,'(a,i3,a)')
     .     '(1x,a,i2,a,1p,',ns1,'(1x,g10.3,'',''))'
          if(minval(dtco(0:ns-1,0:nnreg(0))).ne.1.0_R8 .or.
     &       maxval(dtco(0:ns-1,0:nnreg(0))).ne.1.0_R8) then
            do ireg = 0, nnreg(0)
              cftmp(0:ns-1)=dtco(0:ns-1,ireg)
              call convert_bundle_r
     .         (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
              dtco(0:ns1-1,ireg)=cftmp(0:ns1-1)
              write(99,hlp_frm)
     .         'DTCO(0,',ireg,')=',(dtco(is,ireg),is=0,ns1-1)
            enddo
          endif
          if(minval(dtmo(0:ns-1,0:nnreg(0))).ne.1.0_R8 .or.
     &       maxval(dtmo(0:ns-1,0:nnreg(0))).ne.1.0_R8) then
            do ireg = 0, nnreg(0)
              cftmp(0:ns-1)=dtmo(0:ns-1,ireg)
              call convert_bundle_r
     .         (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
              dtmo(0:ns1-1,ireg)=cftmp(0:ns1-1)
              write(99,hlp_frm)
     .         'DTMO(0,',ireg,')=',(dtmo(is,ireg),is=0,ns1-1)
            enddo
          endif
          write(hlp_frm,'(a,i2,a)')
     .     '(1x,a,1p,',nnreg(0)+1,'(1x,g10.3,'',''))'
          if(minval(dtee(0:nnreg(0))).ne.1.0_R8 .or.
     &       maxval(dtee(0:nnreg(0))).ne.1.0_R8) write(99,hlp_frm)
     .     'DTEE=',(dtee(ireg),ireg=0,nnreg(0))
          if(minval(dtei(0:nnreg(0))).ne.1.0_R8 .or.
     &       maxval(dtei(0:nnreg(0))).ne.1.0_R8) write(99,hlp_frm)
     .     'DTEI=',(dtei(ireg),ireg=0,nnreg(0))
          write(hlp_frm,'(a,i3,a)') '(1x,a,i2,a,',ns1,'l1,a)'
          do ireg = 0, nnreg(0)
            todo = .false.
            do is = 0, ns-1
              todo = todo .or. .not.solveco(is,ireg)
            enddo
            if (todo) then
              call convert_bundle_l(ns1,ns,nsm,is_neutral1,bundle_frac,
     .         solveco(0,ireg))
              write(99,hlp_frm)
     .         'SOLVECO(0,',ireg,')=',(solveco(is,ireg),is=0,ns1-1)
            endif
          enddo
          do ireg = 0, nnreg(0)
            todo = .false.
            do is = 0, ns-1
              todo = todo .or. .not.solvemo(is,ireg)
            enddo
            if (todo) then
              call convert_bundle_l(ns1,ns,nsm,is_neutral1,bundle_frac,
     .         solvemo(0,ireg))
              write(99,hlp_frm)
     .         'SOLVEMO(0,',ireg,')=',(solvemo(is,ireg),is=0,ns1-1)
            endif
          enddo
          write(hlp_frm,'(a,i2,a)') '(1x,a,',nnreg(0)+1,'l1,a)'
          todo = .false.
          do ireg = 0, nnreg(0)
            todo = todo .or. .not.solvepo(ireg)
          enddo
          if (todo) write(99,hlp_frm)
     .      'SOLVEPO=',(solvepo(ireg),is=0,nnreg(0))
          todo = .false.
          do ireg = 0, nnreg(0)
            todo = todo .or. .not.solveee(ireg)
          enddo
          if (todo) write(99,hlp_frm)
     .     'SOLVEEE=',(solveee(ireg),is=0,nnreg(0))
          todo = .false.
          do ireg = 0, nnreg(0)
            todo = todo .or. .not.solveei(ireg)
          enddo
          if (todo) write(99,hlp_frm)
     .     'SOLVEEI=',(solveei(ireg),is=0,nnreg(0))
          if(numerics_time_mod.gt.0.0_R8) then
            write(99,'(1x,a,a,a)')
     .       'NUMERICS_FILENAME=''',trim(numerics_filename),''','
            write(99,'(1x,a,1p,1g10.3,a)')
     .       'NUMERICS_TIME_MOD=',numerics_time_mod,','
            write(99,'(1x,a,1p,1g10.3,a)')
     .       'NUMERICS_TIME_SWITCH=',numerics_time_switch,','
          endif
          if(core_dt_factor.ne.1.0_R8) write(99,'(1x,a,1p,1g10.3,a)')
     .     'CORE_DT_FACTOR=',core_dt_factor,','
          if(core_dt_suppression.ne.1.0_R8) 
     .     write(99,'(1x,a,1p,1g10.3,a)')
     .     'CORE_DT_SUPPRESSION=',core_dt_suppression,','
          if(time_factor_required.ne.0.1_R8)
     .     write(99,'(1x,a,1p,1g10.3,a)')
     .     'TIME_FACTOR_REQUIRED=',time_factor_required,','
          write(99,'(a1)') '/'
          close(99)
        endif
        if (.not.streql(numerics_filename,'b2.numerics.parameters'))
     >   then
          first_pass = .false.
          new_filename =
     .     '../'//numerics_filename(1:lnblnk(numerics_filename))//'.2'
          inquire(file=trim(numerics_filename),exist=todo)
          inquire(file=trim(new_filename),exist=found)
          if (found) write(6,*) trim(new_filename), ' already present!'
        else
          todo = .false.
          found = .true.
        endif
        numerics_namelist = 0
      enddo

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
        write(nout(5),'(a)') 'b2.atomic_physics_rescale.parameters.2'
        write(hlp_frm,'(a,i3,a)') '(1x,a,1p,',ns1,'(1x,g10.3,'',''))'
        open(99,file='../b2.atomic_physics_rescale.parameters.2')
        write(99,'(a23)') '&ATOMIC_PHYSICS_RESCALE'
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
        write(99,'(a1)') '/'
        close(99)
      endif

#ifdef B25_EIRENE
*  ..write new input.dat file if present
      IUNIN = 86
      IUNOT = 87
      if (open_file(IUNIN,'input.dat','')) then
cxpb Now we have read the atomic species data (in FIND_PARAM) and
cxpb start by identifying the correspondence with the new run species
        iss = 0
        DO IATM=1,NATMI
          found = .false.
          is = 0
          do while (.not.found .and. is.lt.ns1)
            found = nint(am1(is)).eq.NMASSA(IATM) .and.
     .              nint(zn1(is)).eq.NCHARA(IATM) .and.
     .              is_neutral1(is)
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
     &           .and. nint(zn1(is)).eq.NCHARA(IATM)
     &           .and. nint(am1(is)).eq.NMASSA(IATM)
              end do
            else
              is = is + 1
            endif
          enddo
          if (.not.found) then
            iatm1(iatm) = 0
            write(*,*) 'Atom species corresponding to ',
     .       trim(TEXTA(IATM)),' not found! Removing...'
          end if
        END DO
        NATMI1 = iss
        if (streql(eirene_format,'iter') .and. same_bundling .and.
     .      nx.eq.nx1 .and. ny.eq.ny1 .and. non_commensurate.eq.0)
     .   goto 8688  ! no need to convert input.dat file
cxpb  Start reading input.dat file
        read(IUNIN,'(A72)',end=8687) ZEILE   ! TXTRUN
        open(unit=IUNOT,file='input.dat.2',status='unknown')
        write(nout(5),'(a)') 'input.dat.2'
        write(*,*) 'Writing input.dat.2'
        write(IUNOT,'(A72)') ZEILE
*** 1. DATA FOR OPERATING MODE
  109   READ (IUNIN,'(A80)') ZEILE
        IF (ZEILE(1:1).EQ.'*') THEN
          WRITE(IUNOT,'(A80)') ZEILE
          GOTO 109
        ELSE
          IF (streql(eirene_format,'old')) THEN
            READ (ZEILE,6666) idum(0:5)  ! NMACH,NMODE,NTIME,NFILE,NITER,lfrm15
            idum(5) = idum(4)
            idum(4) = 1
            idum(6) = 1
            idum(7) = 1
            idum(8) = 0
          ELSE
            READ (ZEILE,6666) idum(0:8)  ! NMACH,NMODE,NTCPU,NFILE,NITER0,NITER,
                                         ! NTIME0,NTIME,NFULLRUN
          END IF
          NMODE = idum(1)
          WRITE(IUNOT,6666) idum(0:8)
        ENDIF
        READ (IUNIN,'(A80)') ZEILE
        IF ((INDEX(ZEILE,'F') + INDEX(ZEILE,'f') + INDEX(ZEILE,'T') +
     .       INDEX(ZEILE,'t')) == 0) THEN
          READ (ZEILE,6666) idum(0:8)  ! IDUM,NOPTM1,NGEOM_USR,NCOUP_INPUT,
                                       ! NSMSTRA,NSTORAM,NGSTAL,NRTAL1,NREAC_ADD
          WRITE(IUNOT,'(A80)') ZEILE
          IF (streql(eirene_format,'facelift')) then
            READ (IUNIN,6665) ldum(0:10) ! NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
                                         ! NLERG,LIDENT,LHABER,NEXVS,NRESCFLUX,nltrimesh
            ldum(14) = ldum(8)
            ldum(15) = ldum(10)
            ldum(7)  = .false.
            ldum(8)  = .false.
            ldum(9)  = .false.
            ldum(10) = .true.
            ldum(11) = .false.
            ldum(12) = .false.
            ldum(13) = .true.
          ELSE
            READ (IUNIN,6665) ldum(0:15) ! NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
                                         ! NLERG,NLIDENT,NLONE,LTSTV,NLDFST,
                                         ! NLOLDRAN,NLCASCAD,NLOCTREE,NLWRMSH,NEXVS,
                                         ! NLTRIMESH
            IF (streql(eirene_format,'old') .or.
     .          streql(eirene_format,'new')) ldum(10) = .true.
          ENDIF
          WRITE(IUNOT,6665) ldum(0:15)
        ELSE
          IF (streql(eirene_format,'facelift')) then
            READ (ZEILE,6665) ldum(0:10) ! NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
                                         ! NLERG,LIDENT,LHABER,NEXVS,NRESCFLUX,nltrimesh
            ldum(14) = ldum(8)
            ldum(15) = ldum(10)
            ldum(7)  = .false.
            ldum(8)  = .false.
            ldum(9)  = .false.
            ldum(10) = .true.
            ldum(11) = .false.
            ldum(12) = .false.
            ldum(13) = .false.
          ELSE
            READ (ZEILE,6665) ldum(0:15) ! NLSCL,NLTEST,NLANA,NLDRFT,NLCRR,
                                         ! NLERG,NLIDENT,NLONE,LTSTV,NLDFST,
                                         ! NLOLDRAN,NLCASCAD,NLOCTREE,NLWRMSH,NEXVS,
                                         ! NLTRIMESH
            IF (streql(eirene_format,'old') .or.
     .          streql(eirene_format,'new')) ldum(10) = .true.
          ENDIF
          WRITE(IUNOT,6665) ldum(0:15)
        ENDIF
        NLWRMSH = ldum(13)
        READ (IUNIN,'(A80)') ZEILE
        I1 = INDEX(ZEILE,'CFILE')
        IF (I1.EQ.0) THEN
          WRITE(IUNOT,'(A80)') ZEILE
          IREAD=1
        ELSE
          DO WHILE (I1 /= 0)
            if (streql(eirene_format,'iter')) then ! no conversion necessary
              WRITE(IUNOT,'(A80)') ZEILE
            else
              I2 = VERIFY(ZEILE(I1+5:),' ') + I1 + 4
              I3 = SCAN(ZEILE(I2+1:),' ')
              HANDLE=REPEAT(' ',6)
              HANDLE(1:I3) = ZEILE(I2:I2+I3-1)
c   cfile card found. Is this one of the permitted external files?
              if (streql(handle,'AMJUEL')) then
                WRITE(IUNOT,'(A47)')
     .           'CFILE AMJUEL Database/AMdata/amjuel.tex        '
              else if (streql(handle,'HYDHEL')) then
                WRITE(IUNOT,'(A47)')
     .           'CFILE HYDHEL Database/AMdata/hydhel.tex        '
              else if (streql(handle,'SPUTER')) then
                WRITE(IUNOT,'(A47)')
     .           'CFILE SPUTER Database/Surfacedata/SPUTER       '
              else if (streql(handle,'TRIM')) then
                WRITE(IUNOT,'(A47)')
     .           'CFILE TRIM   Database/Surfacedata/TRIM/trim.dat'
              else if (streql(handle,'METHAN')) then
                WRITE(IUNOT,'(A47)')
     .           'CFILE METHAN Database/AMdata/methane.tex       '
              else if (streql(handle,'gr_ext')) then
                WRITE(IUNOT,'(A47)')
     .           'CFILE gr_ext Database/graphite_ext.dat         '
              else if (streql(handle,'mo_ext')) then
                WRITE(IUNOT,'(A47)')
     .           'CFILE mo_ext Database/mo_ext.dat               '
              else if (streql(handle,'ADAS')) then
                WRITE(IUNOT,'(A47)')
     .           'CFILE ADAS   modules/adas/adf11                '
              else
                WRITE(IUNOT,'(A80)') ZEILE
              ENDIF
            endif
            READ (IUNIN,'(A80)') ZEILE
            I1 = INDEX(ZEILE,'CFILE')
          END DO
          IREAD=2
        END IF

*** 2. DATA FOR STANDARD MESH
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        ELSE IF (IREAD.EQ.2) THEN
          WRITE(IUNOT,'(A80)') ZEILE
        ENDIF
  201   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 201
        READ (ZEILE,6666) (INDGRD(J),J=1,3)
C  RADIAL MESH
  210   READ (IUNIN,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') THEN
          WRITE(IUNOT,'(A80)') ZEILE
          GOTO 210
        ENDIF
        READ (ZEILE,6665) ldum(0)    ! nlrad
        WRITE(IUNOT,6665) ldum(0)
        IF (ldum(0)) THEN
          READ (IUNIN,6665) ldum(0:7)  ! NLSLB,NLCRC,NLELL,NLTRI,NLPLG,NLFEM,NLTET,NLGEN
          NLTRI = ldum(3)
          NLTET = ldum(6)
          WRITE(IUNOT,6665) ldum(0:7)
          READ (IUNIN,6666) idum(0:4)  ! NR1ST,NRSEP,NRPLG,NPPLG,NRKNOT
          NR1ST = idum(0)
          NPPLG = idum(3)
          if (ldum(0).or.ldum(1).or.ldum(2).or.ldum(3).or.ldum(4)) then
            idum(0) = ny1+1
            IF (INDGRD(1).LE.4) idum(1) = jsep1+1
            idum(2) = nx1 + max(1,nncut*3)
          endif
          idum(5:6) = 0
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
              READ (IUNIN,6666) (NPOINT(1,K),NPOINT(2,K),K=1,idum(3))
              WRITE(IUNOT,6666) (NPOINT(1,K),NPOINT(2,K),K=1,idum(3))
              DO 212 I=1,NR1ST
                READ (IUNIN,6664) (XPOL(I,J),YPOL(I,J),J=1,idum(2))
                WRITE(IUNOT,6664) (XPOL(I,J),YPOL(I,J),J=1,idum(2))
212           CONTINUE
              deallocate(NPOINT)
              deallocate(XPOL,YPOL)
            ENDIF
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
          ELSEIF (INDGRD(1).EQ.6) THEN
            READ (IUNIN,'(A80)') ZEILE
            WRITE(IUNOT,'(A80)') ZEILE
            IPOS1=INDEX(ZEILE,'T')
            IPOS2=INDEX(ZEILE,'F')
            IF (IPOS1.GT.0.OR.IPOS2.GT.0) GOTO 222
          ENDIF
        ENDIF
C  POLOIDAL MESH
  220   READ (IUNIN,'(A80)') ZEILE  ! NLPOL
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 220
  222   READ (IUNIN,6665) ldum(0:2) ! NLPLY,NLPLA,NLPLP
        WRITE(IUNOT,6665) ldum(0:2)
        READ (IUNIN,6666) idum(0:3) ! NP2ND,NPSEP,NPPLA,NPPER
        NP2ND = idum(0)
        idum(0) = nx1 + max(1,nncut*3)
        WRITE(IUNOT,6666) idum(0:3)
        IF (INDGRD(2).LE.5) THEN
          READ (IUNIN,'(A72)') ZEILE ! YIA,YGA,YAA,YYA
          WRITE(IUNOT,'(A72)') ZEILE
        ELSEIF (INDGRD(2).EQ.6) THEN
        ENDIF
C  TOROIDAL MESH
  230   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 230
        READ (IUNIN,6665) ldum(0:2) ! NLTRZ,NLTRA,NLTRT
        WRITE(IUNOT,6665) ldum(0:2)
        READ (IUNIN,6666) idum(0:3) ! NT3RD,NTSEP,NTTRA,NTPER
        NT3RD = idum(0)
        WRITE(IUNOT,6666) idum(0:3)
        IF (INDGRD(3).LE.5) THEN
          READ (IUNIN,'(A72)') ZEILE ! ZIA,ZGA,ZAA,ZZA,ROA
          WRITE(IUNOT,'(A72)') ZEILE
        ELSEIF (INDGRD(3).EQ.6) THEN
        ENDIF
C  MESH MULTIPLICATION
  240   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 240
        IREAD=1
        READ (ZEILE,6665) ldum(0)
        IREAD=0
        IF (ldum(0)) THEN
          READ (IUNIN,6666) idum(0)
          WRITE(IUNOT,6666) idum(0)
          allocate(VOLCOR(idum(0)))
          READ (IUNIN,6664) (VOLCOR(I),I=1,idum(0))
          WRITE(IUNOT,6664) (VOLCOR(I),I=1,idum(0))
          deallocate(VOLCOR)
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        ELSE
  241     READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
          IPOS0=INDEX(ZEILE,'*')
          IPOS1=INDEX(ZEILE,'T')
          IPOS2=INDEX(ZEILE,'F')
          IF (IPOS0.EQ.0.AND.IPOS1.EQ.0.AND.IPOS2.EQ.0) GOTO 241
        ENDIF
        IREAD=1
C  ADDITIONAL CELLS OUTSIDE STANDARD MESH
  250   IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        ENDIF
        IF (ZEILE(1:1) .EQ. '*') THEN
          IREAD=0
          GOTO 250
        ENDIF
        READ (ZEILE,6665) ldum(0)
        IREAD=0
        IF (ldum(0)) THEN
          READ (IUNIN,6666) idum(0) ! NRADD
          WRITE(IUNOT,6666) idum(0)
          allocate(VOLADD(idum(0)))
          READ (IUNIN,6664) (VOLADD(I),I=1,idum(0))
          WRITE(IUNOT,6664) (VOLADD(I),I=1,idum(0))
          deallocate(VOLADD)
        ELSE
C  FIND START OF NEXT INPUT BLOCK: 3A
  252     READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
          IF (ZEILE(1:3) .NE. '***') GOTO 252
          IREAD=1
        ENDIF
*** 3A. DATA FOR NON DEFAULT STANDARD SURFACES
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        ENDIF
  310   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 310
        IREAD=1
        READ (ZEILE,6666) NSTSI
        IREAD=0
        allocate(ILIIN(NLIM+NSTS))
        DO 311 ISTS=1,NSTSI
          IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A80)') ZEILE ! TXTSFL(NLJ)
            WRITE(IUNOT,'(A80)') ZEILE
          ELSE
            IREAD=0
          ENDIF
          READ (IUNIN,6666) idum(0:8) ! JDUMMY,IDIMP,INUMP(ISTS,IDIMP),IRPTA1,
                                      ! IRPTE1,IRPTA2,IRPTE2,IRPTA3,IRPTE3
          IDIMP  = idum(1)
          INUMP  = idum(2)
          IRPTA1 = idum(3)
          IRPTE1 = idum(4)
          IRPTA2 = idum(5)
          IRPTE2 = idum(6)
          IRPTA3 = idum(7)
          IRPTE3 = idum(8)
C  OLD INPUT VERSION BEGIN
          IF (IDIMP.EQ.1.AND.IRPTA1.NE.IRPTE1) THEN
            IRPTA2=IRPTA1
            IRPTE2=IRPTE1
            IRPTA1=INUMP
            IRPTE1=INUMP
          ELSE IF (IDIMP.EQ.1.AND.IRPTA1.EQ.IRPTE1.AND.
     .             IRPTA2+IRPTE2+IRPTA3+IRPTE3.EQ.0) THEN
            IRPTA2=IRPTA1
            IRPTE2=IRPTE1
            IRPTA1=INUMP
            IRPTE1=INUMP
          ELSEIF (IDIMP.EQ.1.AND.IRPTA1.EQ.0.AND.IRPTE1.EQ.0) THEN
            IRPTA1=INUMP
            IRPTE1=INUMP
          ELSEIF (IDIMP.EQ.2.AND.IRPTA2.EQ.0.AND.IRPTE2.EQ.0) THEN
            IRPTA2=INUMP
            IRPTE2=INUMP
          ELSEIF (IDIMP.EQ.3.AND.IRPTA3.EQ.0.AND.IRPTE3.EQ.0) THEN
            IRPTA3=INUMP
            IRPTE3=INUMP
          ENDIF
          IF (IDIMP.EQ.1) THEN
            idum(2) = convert_eirene(non_commensurate,
     .       1,INUMP,(IRPTA2+IRPTE2)/2,nx,ny,nx1,ny1,4)
            idum(3) = idum(2)
            idum(4) = idum(2)
            idum(5) = convert_eirene(non_commensurate,
     .       2,INUMP,IRPTA2,nx,ny,nx1,ny1,1)
            idum(6) = convert_eirene(non_commensurate,
     .       2,INUMP,IRPTE2,nx,ny,nx1,ny1,3)
          ELSE IF (IDIMP.EQ.2) THEN
            idum(2) = convert_eirene(non_commensurate,
     .       2,(IRPTA1+IRPTE1)/2,INUMP,nx,ny,nx1,ny1,4)
            idum(3) = convert_eirene(non_commensurate,
     .       1,IRPTA1,INUMP,nx,ny,nx1,ny1,0)
            idum(4) = convert_eirene(non_commensurate,
     .       1,IRPTE1,INUMP,nx,ny,nx1,ny1,2)
            idum(5) = idum(2)
            idum(6) = idum(2)
          ELSE IF (IDIMP.EQ.3) THEN
            idum(3) = convert_eirene(non_commensurate,
     .       1,IRPTA1,IRPTA2,nx,ny,nx1,ny1,0)
            idum(4) = convert_eirene(non_commensurate,
     .       1,IRPTE1,IRPTE2,nx,ny,nx1,ny1,2)
            idum(5) = convert_eirene(non_commensurate,
     .       2,IRPTA1,IRPTA2,nx,ny,nx1,ny1,1)
            idum(6) = convert_eirene(non_commensurate,
     .       2,IRPTE1,IRPTE2,nx,ny,nx1,ny1,3)
            idum(7) = IRPTA3
            idum(8) = IRPTE3
          ENDIF
          WRITE(IUNOT,6666) idum(0:8)
          READ (IUNIN,6666) idum(0:9) ! ILIIN(NLJ),ILSIDE(NLJ),ILSWCH(NLJ),
                                      ! ILEQUI(NLJ),ILTOR(NLJ),ILCOL(NLJ),
                                      ! ILFIT(NLJ),ILCELL(NLJ),ILBOX(NLJ),
                                      ! ILPLG(NLJ)
          ILIIN(ISTS) = idum(0)
          ILEQUI = idum(3)
          IF (IDIMP.EQ.1.AND.ILEQUI.NE.0) THEN
            idum(3) = convert_eirene(non_commensurate,1,
     .       ILEQUI,(IRPTA1+IRPTE1)/2,nx,ny,nx1,ny1,4)
          ELSE IF (IDIMP.EQ.2.AND.ILEQUI.NE.0) THEN
            idum(3) = convert_eirene(non_commensurate,2,
     .       (IRPTA1+IRPTE1)/2,ILEQUI,nx,ny,nx1,ny1,4)
          ENDIF
          WRITE(IUNOT,6666) idum(0:9)
  312     READ (IUNIN,'(A80)') ZEILE
          IREAD=1
          IF (ZEILE(1:1).EQ.'*') THEN
            GOTO 314
          ELSEIF (ZEILE(1:8).EQ.'SURFMOD_') THEN
            IREAD=0
            GOTO 314
          ELSE
            READ (ZEILE,6666) idum(0:4) ! ILREF(NLJ),ILSPT(NLJ),ISRS(1,NLJ),
                                        ! ISRC(1,NLJ),LCHSPNWL(1,NLJ)
            IF (idum(2).gt.0 .and. idum(2).le.NATMI)
     >        idum(2) = iatm1(idum(2))
            IF (idum(3).gt.0 .and. idum(3).le.NATMI) THEN
              idum(3) = iatm1(idum(3))
            ELSE IF (idum(3).gt.NATMI .and. idum(3).le.NATMI+NMOLI) THEN
              idum(3) = idum(3) - NATMI + NATMI1
            ENDIF
            IF (streql(eirene_format,'juelich') .or.
     .          streql(eirene_format,'old') .or.
     .          streql(eirene_format,'new')) idum(4) = 1
            WRITE(IUNOT,6666) idum(0:4)
            IREAD=0
            IF (streql(eirene_format,'old')) then
              READ (IUNIN,6664) rdum(0:5) ! ZNML(NLJ),EWALL(NLJ),EWBIN(NLJ),
                                          ! TRANSP(1,NLJ),TRANSP(2,NLJ),RECYCC(NLJ)
              WRITE(IUNOT,6664) rdum(0:4),0.0_R8
              READ (IUNIN,6664) rdum(6:11) ! RECYCF(1,NLJ),RECYCT(1,NLJ),RECYCS(1,NLJ),
                                           ! EXPPL(1,NLJ),EXPEL(1,NLJ),EXPIL(1,NLJ)
              WRITE(IUNOT,6664) rdum(6:7),0.0_R8,rdum(9:11)
              WRITE(IUNOT,6664) rdum(8),rdum(5),1.0_R8,0.0_R8,0.0_R8
            ELSE
              READ (IUNIN,'(A72)') ZEILE ! ZNML(NLJ),EWALL(NLJ),EWBIN(NLJ),
                                         ! TRANSP(1,NLJ),TRANSP(2,NLJ),FSHEAT(NLJ)
              WRITE(IUNOT,'(A72)') ZEILE
              READ (IUNIN,'(A72)') ZEILE ! RECYCF(1,NLJ),RECYCT(1,NLJ),RECPRM(1,NLJ),
                                         ! EXPPL(1,NLJ),EXPEL(1,NLJ),EXPIL(1,NLJ)
              WRITE(IUNOT,'(A72)') ZEILE
            ENDIF
            READ (IUNIN,'(A80)') ZEILE
            IREAD=1
C  READ ONE MORE LINE FOR NON-DEFAULT SPUTTER MODEL
            IF (ZEILE(1:1).NE.'*') THEN
              IF (streql(eirene_format,'facelift')) THEN
                READ (ZEILE,6664) rdum(0:1), ! RECYCS(1,NLJ),RECYCC(1,NLJ),
     .                     rdum(4), rdum(5), ! ECHSPTG(1,NLJ),FCHSPTG(1,NLJ),
     .                     rdum(2)           ! SPTPRM(1,NLJ)
                rdum(3) = 0.0_R8
              ELSE IF (streql(eirene_format,'iter')) THEN
                READ (ZEILE,6664) rdum(0:4)  ! RECYCS(1,NLJ),RECYCC(1,NLJ),
                                             ! SPTPRM(1,NLJ),ESPUTS(1,NLJ),
                                             ! ESPUTC(1,NLJ)
              ELSE
                READ (ZEILE,6664) rdum(0:2)  ! RECYCS(1,NLJ),RECYCC(1,NLJ),
                                             ! SPTPRM(1,NLJ)
                rdum(3) = 0.0_R8
                rdum(4) = 0.0_R8
              ENDIF
              IREAD=0
              WRITE(IUNOT,6664) rdum(0:4)
              GOTO 311
            ENDIF
          ENDIF
  314     CONTINUE
          WRITE(IUNOT,'(A80)') ZEILE
  311   CONTINUE
*** 3B. DATA FOR ADDITIONAL SURFACES
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        ENDIF
  350   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:1).EQ.'*') GOTO 350
        IREAD=0
        READ (ZEILE,6666) NLIMI
        IF (NLIMI.GT.0) THEN
          allocate(IHELP(NLIMI),IGJUM0(NLIMI))
          IGJUM0(1:NLIMI) = 0
          IHELP(1:NLIMI) = 0
  351     READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
          IREAD=1
          IF (ZEILE(1:3).EQ.'CH0') THEN
            IREAD=0
            CALL EIRENE_DEKEY (ZEILE(4:80),IHELP,1,1,1,NLIMI)
            GOTO 351
          ENDIF
          DO 352 I=1,NLIMI
            IGJUM0(I)=IHELP(I)
  352     CONTINUE
          allocate(RLB(NLIMI))
        ENDIF
        DO 360 I=1,NLIMI
          IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A80)') ZEILE ! TXTSFL(I)
            WRITE(IUNOT,'(A80)') ZEILE
          ENDIF
          IREAD=0
          IF (IGJUM0(I).NE.0) THEN
  361       READ (IUNIN,'(A80)') ZEILE
            WRITE(IUNOT,'(A80)') ZEILE
            IREAD=1
            IF (ZEILE(1:1).EQ.'*') GOTO 369
            IF (ZEILE(1:9).EQ.'TRANSFORM') GOTO 368
            GOTO 361
          ENDIF
  362     READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
          IF (ZEILE(1:3).EQ.'CH1') THEN
            GOTO 362
          ELSEIF (ZEILE(1:3).EQ.'CH2') THEN
            GOTO 362
          ELSE
            READ (ZEILE,6664) RLB(I),rdum(1:3) ! SAREA(I),RLWMN(I),RLWMX(I)
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
  367     READ (IUNIN,'(A80)') ZEILE
          IREAD=1
          IF (ZEILE(1:1).EQ.'*') THEN
            WRITE(IUNOT,'(A80)') ZEILE
            GOTO 369
          ELSEIF (ZEILE(1:9).EQ.'TRANSFORM') THEN
            WRITE(IUNOT,'(A80)') ZEILE
            GOTO 368
          ELSEIF (ZEILE(1:8).EQ.'SURFMOD_') THEN
            WRITE(IUNOT,'(A80)') ZEILE
            IREAD=0
            GOTO 368
          ELSE
            READ (ZEILE,6666) idum(0:4) ! ILREF(I),ILSPT(I),ISRS(1,I),ISRC(1,I),
                                        ! LCHSPNWL(1,I)
            IF (idum(2).gt.0 .and. idum(2).le.NATMI)
     >       idum(2) = iatm1(idum(2))
            IF (idum(3).gt.0 .and. idum(3).le.NATMI)
     >       idum(3) = iatm1(idum(3))
            IF (streql(eirene_format,'juelich') .or.
     .          streql(eirene_format,'old') .or.
     .          streql(eirene_format,'new')) idum(4) = 1
            WRITE(IUNOT,6666) idum(0:4)
            IREAD=0
            IF (streql(eirene_format,'old')) then
              READ (IUNIN,6664) rdum(0:5) ! ZNML(NLJ),EWALL(NLJ),EWBIN(NLJ),
                                          ! TRANSP(1,NLJ),TRANSP(2,NLJ),RECYCC(NLJ)
              WRITE(IUNOT,6664) rdum(0:4),0.0_R8
              READ (IUNIN,6664) rdum(6:11) ! RECYCF(1,NLJ),RECYCT(1,NLJ),RECYCS(1,NLJ),
                                           ! EXPPL(1,NLJ),EXPEL(1,NLJ),EXPIL(1,NLJ)
              WRITE(IUNOT,6664) rdum(6:7),0.0_R8,rdum(9:11)
              WRITE(IUNOT,6664) rdum(8),rdum(5),1.0_R8,0.0_R8,0.0_R8
            ELSE
              READ (IUNIN,'(A72)') ZEILE ! ZNML(I),EWALL(I),EWBIN(I),
                                         ! TRANSP(1,I),TRANSP(2,I),FSHEAT(I)
              WRITE(IUNOT,'(A72)') ZEILE
              READ (IUNIN,'(A72)') ZEILE ! RECYCF(1,I),RECYCT(1,I),RECPRM(1,I),
                                         ! EXPPL(1,I),EXPEL(1,I),EXPIL(1,I)
              WRITE(IUNOT,'(A72)') ZEILE
            ENDIF
            READ (IUNIN,'(A80)') ZEILE
            IREAD=1
            IF (ZEILE(1:1).NE.'*'.AND.ZEILE(1:9).NE.'TRANSFORM') THEN
              IF (streql(eirene_format,'facelift')) THEN
                READ (ZEILE,6664) rdum(0:1), ! RECYCS(1,NLJ),RECYCC(1,NLJ),
     .                     rdum(4), rdum(5), ! ECHSPTG(1,NLJ),FCHSPTG(1,NLJ),
     .                     rdum(2)           ! SPTPRM(1,NLJ)
                rdum(3) = 0.0_R8
              ELSE IF (streql(eirene_format,'iter')) THEN
                READ (ZEILE,6664) rdum(0:4)  ! RECYCS(1,NLJ),RECYCC(1,NLJ),
                                             ! SPTPRM(1,NLJ),ESPUTS(1,NLJ),
                                             ! ESPUTC(1,NLJ)
              ELSE
                READ (ZEILE,6664) rdum(0:2)  ! RECYCS(1,NLJ),RECYCC(1,NLJ),
                                             ! SPTPRM(1,NLJ)
                rdum(3) = 0.0_R8
                rdum(4) = 0.0_R8
              ENDIF
              IREAD=0
              WRITE(IUNOT,6664) rdum(0:4)
            ELSE
              WRITE(IUNOT,'(A80)') ZEILE
            ENDIF
          ENDIF
  368     IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A80)') ZEILE
            WRITE(IUNOT,'(A80)') ZEILE
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
        IF (NLIMI.GT.0) THEN
          deallocate(RLB)
          deallocate(ILIIN,IHELP,IGJUM0)
        ENDIF
C  READ DATA FOR SPECIES SPECIFICATION AND ATOMIC PHYSICS MODULE
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A72)') ZEILE
          WRITE(IUNOT,'(A72)') ZEILE
        ENDIF
        IREAD=0
        READ (IUNIN,'(A72)') ZEILE
        WRITE(IUNOT,'(A72)') ZEILE
        READ (IUNIN,*) idum(0)  ! NREACI
        allocate(IREAC1(idum(0)))
        allocate(REAC_LINE(3*idum(0)+1))
        IREAC1 = -1
        NREAC1 = 0
        ILINE = 0
  411   READ (IUNIN,'(A80)') ZEILE
        IF (ZEILE(1:1).NE.'*') THEN
          IF (INDEX(ZEILE,'PHOTON').EQ.0) THEN
            READ (ZEILE,66661)
     .            I,FILNAM,H123,REAC,CRC,idum(0:1),rdum(0),RMN,RMX
!  CORRECT FOR OLD REACTIONS NO LONGER SUPPORTED
!  Remove AMJUEL H.6 reactions
!  Remove AMJUEL H.9 0.2/0.3/0.4 reactions
            skip_reac = .false.
            IF(streql(FILNAM,'AMJUEL')) THEN
              IF (streql(REAC,'0.1')) skip_reac = .true.
              IF (streql(H123,'H.6')) skip_reac = .true.
              IF (streql(H123,'H.9')) THEN
                IF (streql(REAC,'0.2')) skip_reac = .true.
                IF (streql(REAC,'0.3')) skip_reac = .true.
                IF (streql(REAC,'0.4')) skip_reac = .true.
              END IF
            END IF
            IF (.not.skip_reac) THEN
              IF (IREAC1(I).EQ.-1) THEN ! To account for reaction with multiple rates
                NREAC1=NREAC1+1
                IREAC1(I)=NREAC1
              END IF
              ILINE=ILINE+1
              IF (streql(FILNAM,'AMJUEL')) THEN
                IF (streql(REAC,'0.2')) REAC='0.2T'
                IF (streql(REAC,'0.3')) REAC='0.3T'
                IF (streql(REAC,'0.4')) REAC='0.4T'
                IF (streql(REAC,'0.5')) REAC='0.5T'
                IF (streql(REAC,'0.6')) REAC='0.6T'
                IF (streql(REAC,'0.7')) REAC='0.7T'
                IF (streql(REAC,'0.8')) REAC='0.8T'
                IF (streql(REAC,'7.0c') .and. streql(H123,'H.11'))
     .           REAC='7.0a'
              END IF
              WRITE(REAC_LINE(ILINE),66661) IREAC1(I),
     .            FILNAM,H123,REAC,CRC,idum(0:1),rdum(0),RMN,RMX
            END IF
          ELSE
            NREAC1 = NREAC1 + 1
            ZEILE=ADJUSTL(ZEILE)
            I1=SCAN(ZEILE,' ')
            READ (ZEILE(1:I1-1),*) I
!  READ FILNAM
            IANF=I1+VERIFY(ZEILE(I1+1:),' ')
            IEND=IANF+SCAN(ZEILE(IANF+1:),' ')-1
            FILNAM(1:IEND-IANF+1)=ZEILE(IANF:IEND)
!  READ H123
            IANF=IEND+VERIFY(ZEILE(IEND+1:),' ')
            IEND=IANF+SCAN(ZEILE(IANF+1:),' ')-1
            H123='    '
            H123(1:IEND-IANF+1)=ZEILE(IANF:IEND)
!  READ REAC2
            IANF=IEND+VERIFY(ZEILE(IEND+1:),' ')
            IEND=IANF+SCAN(ZEILE(IANF+1:),' ')-1
            REAC2=REPEAT(' ',50)
            REAC2(1:IEND-IANF+1)=ZEILE(IANF:IEND)
            REACSTR=REPEAT(' ',11)
            IANF=VERIFY(REAC2,' ')
            IANF=MAX(1,IANF) !VK
            I=INDEX(REAC2(IANF:),' ')-1
            IF (I.LT.0) I=lnblnk(REAC2(IANF:))
            REACSTR(2:I+1)=REAC2(IANF:I+IANF-1)
!  READ CRC
            IANF=IEND+VERIFY(ZEILE(IEND+1:),' ')
            IEND=IANF+SCAN(ZEILE(IANF+1:),' ')-1
            CRC(1:IEND-IANF+1)=ZEILE(IANF:IEND)
            READ (ZEILE(IEND+1:),*) idum(0:1),rdum(0),RMN,RMX
            ILINE = ILINE+1
            WRITE(REAC_LINE(ILINE),'(A80)') ZEILE
          ENDIF
          IF (INDEX(ZEILE,'ADAS') .NE. 0) THEN
            READ (IUNIN,'(4X,A2,1X,I3)') ELNAME,idum(0)
            ILINE = ILINE+1
            WRITE(REAC_LINE(ILINE),'(4X,A2,1X,I3)') ELNAME,idum(0)
          ENDIF
          IF (INDEX(H123,'P.').EQ.0) THEN  !PVK
            IF (RMN.GT.0.0_R8) THEN
              READ (IUNIN,66664) idum(0),(rdum(I),I=1,3) ! IFEXMN(IR,J),(FPARM(IR,I,J),I=1,3)
              if (.not.skip_reac) then
                ILINE = ILINE+1
                WRITE(REAC_LINE(ILINE),66664) idum(0),(rdum(I),I=1,3)
              end if
            ENDIF
            IF (RMX.GT.0.0_R8) THEN
              READ (IUNIN,66664) idum(0),(rdum(I),I=4,6) ! IFEXMN(IR,J),(FPARM(IR,I,J),I=4,6)
              if (.not.skip_reac) then
                ILINE = ILINE+1
                WRITE(REAC_LINE(ILINE),66664) idum(0),(rdum(I),I=4,6)
              end if
            ENDIF
          ENDIF
          IF (INDEX(ZEILE,'CONST') .NE. 0) THEN
            I=INDEX(ZEILE(16:),'FT')
            IF (I /= 0) THEN
              READ (REACSTR(I+2:),*) idum(0) ! IFTFLG(IR,IFLG)
            END IF
            IF (MOD(idum(0),100) == 10) THEN
              ILINE = ILINE+1
              READ (IUNIN,6664) rdum(0) ! CREACD(1,1)
              WRITE(REAC_LINE(ILINE),6664) rdum(0)
            ELSE
              ILINE = ILINE+1
              READ (IUNIN,6664) (rdum(I),I=1,9) ! (CREACD(IC,1),IC=1,9)
              WRITE(REAC_LINE(ILINE),6664) (rdum(I),I=1,6)
              ILINE = ILINE+1
              WRITE(REAC_LINE(ILINE),6664) (rdum(I),I=7,9)
            ENDIF
          ENDIF
          GOTO 411
        ELSE
          IF (ILINE.EQ.0) THEN ! Block 4 header lines
            WRITE(IUNOT,'(A80)') ZEILE
          ELSE ! Block 4a header lines
            ILINE = ILINE+1
            WRITE(REAC_LINE(ILINE),'(A80)') ZEILE
          ENDIF
        ENDIF
! We now write the new list of reactions
        WRITE(IUNOT,6666) NREAC1
        DO I = 1, ILINE
          WRITE(IUNOT,'(A80)') REAC_LINE(I)
        END DO
        deallocate(reac_line)
*4A.   NEUTRAL ATOMS SPECIES CARDS, NATMI SPECIES
Cxpb We must read up to section 5 first, to identify ion species
Cxpb and bundles, compute new bulk ion species indices and only
Cxpb then rewrite the entire block
!! FIXME : Need to consider cases where a molecule or a test ion disappear
!! FIXME : because their constituents atoms are removed
  420   CONTINUE
        READ (IUNIN,*) NATMI   ! last read before possible reassignments
        allocate(NPRTA(NATMI),
     .           ISRFA(NATMI),ISRTA(NATMI),ISRSA(NATMI),
     .           NRCA(NATMI))
        allocate(NFOLA(NATMI),NGENA(NATMI),NHSTA(NATMI))
        DO 421 IATM=1,NATMI
          READ (IUNIN,66666) I,TEXTA(IATM),NMASSA(IATM),
     .                         NCHARA(IATM),NPRTA(IATM),idum(3),
     .                         ISRFA(IATM),ISRTA(IATM),ISRSA(IATM),
     .                         NRCA(IATM),NFOLA(IATM),NGENA(IATM),
     .                         NHSTA(IATM)
          DO 422 K=1,NRCA(IATM)
            IF (ISRSA(IATM).LT.3) THEN
              READ (IUNIN,6666) idum(0:6)
            ELSE IF (ISRSA(IATM).EQ.3) THEN
              READ (IUNIN,6666) idum(0:7)
            ELSE IF (ISRSA(IATM).EQ.4) THEN
              READ (IUNIN,6666) idum(0:8)
            ENDIF
            READ (IUNIN,6664) rdum(0:5)
422       CONTINUE
421     CONTINUE
        DO 4211 IATM=1,NATMI
          backspace(IUNIN)
          DO 4221 K=1,NRCA(IATM)
            backspace(IUNIN)
            backspace(IUNIN)
 4221     CONTINUE
 4211   CONTINUE
        NNREAC=max(1,maxval(NRCA(1:NATMI)))
        allocate(IREACA(NATMI,NNREAC),IBULKA(NATMI,NNREAC),
     .           ISCD1A(NATMI,NNREAC),ISCD2A(NATMI,NNREAC),
     .           ISCD3A(NATMI,NNREAC),ISCD4A(NATMI,NNREAC),
     .           ISCDEA(NATMI,NNREAC))
        allocate(IESTMA(NATMI,NNREAC),IBGKA(NATMI,NNREAC))
        allocate(EELECA(NATMI,NNREAC),EBULKA(NATMI,NNREAC),
     .           ESCD1A(NATMI,NNREAC),ESCD2A(NATMI,NNREAC),
     .           FREACA(NATMI,NNREAC),FLDLMA(NATMI,NNREAC))
        DO 4212 IATM=1,NATMI
          READ (IUNIN,'(A72)') ZEILE
          DO 4222 K=1,NRCA(IATM)
            IF (ISRSA(IATM).LT.3) THEN
              READ (IUNIN,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                          ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                          ISCDEA(IATM,K),IESTMA(IATM,K),
     .                          IBGKA(IATM,K)
            ELSE IF (ISRSA(IATM).EQ.3) THEN
              READ (IUNIN,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                          ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                          ISCD3A(IATM,K),
     .                          ISCDEA(IATM,K),IESTMA(IATM,K),
     .                          IBGKA(IATM,K)
            ELSE IF (ISRSA(IATM).EQ.4) THEN
              READ (IUNIN,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                          ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                          ISCD3A(IATM,K),ISCD4A(IATM,K),
     .                          ISCDEA(IATM,K),IESTMA(IATM,K),
     .                          IBGKA(IATM,K)
            ENDIF
            READ (IUNIN,6664) EELECA(IATM,K),EBULKA(IATM,K),
     .                        ESCD1A(IATM,K),ESCD2A(IATM,K),
     .                        FREACA(IATM,K),FLDLMA(IATM,K)
 4222     CONTINUE
          J=0
          DO K=1,NRCA(IATM)
            IF (IREAC1(IREACA(IATM,K)).GE.1) THEN
              J = J + 1
              IREACA(IATM,J) = IREAC1(IREACA(IATM,K))
              IBULKA(IATM,J) = IBULKA(IATM,K)
              ISCD1A(IATM,J) = ISCD1A(IATM,K)
              ISCD2A(IATM,J) = ISCD2A(IATM,K)
              ISCD3A(IATM,J) = ISCD3A(IATM,K)
              ISCD4A(IATM,J) = ISCD4A(IATM,K)
              ISCDEA(IATM,J) = ISCDEA(IATM,K)
              IESTMA(IATM,J) = IESTMA(IATM,K)
              IBGKA(IATM,J) = IBGKA(IATM,K)
              IF (B2IDEZ(ISCDEA(IATM,J),5,5).eq.3) THEN
                IR = nint(EELECA(IATM,K))
                EELECA(IATM,J) = dble(IREAC1(IR))
              ELSE
                EELECA(IATM,J) = EELECA(IATM,K)
              END IF
              EBULKA(IATM,J) = EBULKA(IATM,K)
              ESCD1A(IATM,J) = ESCD1A(IATM,K)
              ESCD2A(IATM,J) = ESCD2A(IATM,K)
              FREACA(IATM,J) = FREACA(IATM,K)
              FLDLMA(IATM,J) = FLDLMA(IATM,K)
            END IF
          END DO
          NRCA(IATM) = J
 4212   CONTINUE
*4B.   NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES
        READ (IUNIN,'(A72)') TXTMOL
  430   READ (IUNIN,*) NMOLI
        allocate(TEXTM(NMOLI))
        allocate(NMASSM(NMOLI),NCHARM(NMOLI),NPRTM(NMOLI),
     .           ISRFM(NMOLI),ISRTM(NMOLI),ISRSM(NMOLI),
     .           NRCM(NMOLI))
        allocate(NFOLM(NMOLI),NGENM(NMOLI),NHSTM(NMOLI))
        DO 431 IMOL=1,NMOLI
          IF (streql(eirene_format,'juelich') .or.
     .        streql(eirene_format,'iter')) THEN
            READ (IUNIN,66666) I,TEXTM(IMOL),NMASSM(IMOL),NCHARM(IMOL),
     .                         NPRTM(IMOL),idum(3),
     .                         ISRFM(IMOL),ISRTM(IMOL),ISRSM(IMOL),
     .                         NRCM(IMOL),NFOLM(IMOL),NGENM(IMOL),
     .                         NHSTM(IMOL),LKINDM(IMOL)
          ELSE IF (streql(eirene_format,'facelift') .or.
     .             streql(eirene_format,'new')) THEN
            READ (IUNIN,66666) I,TEXTM(IMOL),NMASSM(IMOL),NCHARM(IMOL),
     .                         NPRTM(IMOL),idum(3),
     .                         ISRFM(IMOL),ISRTM(IMOL),ISRSM(IMOL),
     .                         NRCM(IMOL),NFOLM(IMOL),NGENM(IMOL)
            NHSTM(IMOL) = 0
            LKINDM(IMOL) = LMOLSCL(IMOL)
          ELSE IF (streql(eirene_format,'old')) THEN
            READ (IUNIN,66666) I,TEXTM(IMOL),NMASSM(IMOL),NCHARM(IMOL),
     .                         NPRTM(IMOL),idum(3),
     .                         ISRFM(IMOL),ISRTM(IMOL),ISRSM(IMOL),
     .                         NRCM(IMOL),LKINDM(IMOL)
            NFOLM(IMOL) = 0
            NGENM(IMOL) = 0
            NHSTM(IMOL) = 0
          ENDIF
          IF (ISRTM(IMOL).LT.0) ISRTM(IMOL) = NMOLI+1
          DO 432 K=1,NRCM(IMOL)
            IF (ISRSM(IMOL).LT.3) THEN
              READ (IUNIN,6666) idum(0:6) ! IREACM(IMOL,K),IBULKM(IMOL,K),
                                          ! ISCD1M(IMOL,K),ISCD2M(IMOL,K),
                                          ! ISCDEM(IMOL,K),IESTMM(IMOL,K),
                                          ! IBGKM(IMOL,K)
            ELSE IF (ISRSM(IMOL).EQ.3) THEN
              READ (IUNIN,6666) idum(0:7) ! IREACM(IMOL,K),IBULKM(IMOL,K),
                                          ! ISCD1M(IMOL,K),ISCD2M(IMOL,K),
                                          ! ISCD3M(IMOL,K),
                                          ! ISCDEM(IMOL,K),IESTMM(IMOL,K),
                                          ! IBGKM(IMOL,K)
            ELSE IF (ISRSM(IMOL).EQ.4) THEN
              READ (IUNIN,6666) idum(0:8) ! IREACM(IMOL,K),IBULKM(IMOL,K),
                                          ! ISCD1M(IMOL,K),ISCD2M(IMOL,K),
                                          ! ISCD3M(IMOL,K),ISCD4M(IMOL,K),
                                          ! ISCDEM(IMOL,K),IESTMM(IMOL,K),
                                          ! IBGKM(IMOL,K)
            END IF
            READ (IUNIN,6664) rdum(0:4) ! EELECM(IMOL,K),EBULKM(IMOL,K),
                                        ! ESCD1M(IMOL,K),ESCD2M(IMOL,K),
                                        ! FREACM(IMOL,K)
432       CONTINUE
431     CONTINUE
        DO 4311 IMOL=1,NMOLI
          backspace(IUNIN)
          DO 4321 K=1,NRCM(IMOL)
            backspace(IUNIN)
            backspace(IUNIN)
 4321     CONTINUE
 4311   CONTINUE
        NNREAC=max(1,maxval(NRCM(1:NMOLI)))
        allocate(IREACM(NMOLI,NNREAC),IBULKM(NMOLI,NNREAC),
     .           ISCD1M(NMOLI,NNREAC),ISCD2M(NMOLI,NNREAC),
     .           ISCD3M(NMOLI,NNREAC),ISCD4M(NMOLI,NNREAC),
     .           ISCDEM(NMOLI,NNREAC))
        allocate(IESTMM(NMOLI,NNREAC),IBGKM(NMOLI,NNREAC))
        allocate(EELECM(NMOLI,NNREAC),EBULKM(NMOLI,NNREAC),
     .           ESCD1M(NMOLI,NNREAC),ESCD2M(NMOLI,NNREAC),
     .           FREACM(NMOLI,NNREAC))
        DO 4312 IMOL=1,NMOLI
          READ (IUNIN,'(A72)') ZEILE
          DO 4322 K=1,NRCM(IMOL)
            IF (ISRSM(IMOL).NE.3) THEN
              READ (IUNIN,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                          ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                          ISCDEM(IMOL,K),IESTMM(IMOL,K),
     .                          IBGKM(IMOL,K)
            ELSE IF (ISRSM(IMOL).EQ.3) THEN
              READ (IUNIN,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                          ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                          ISCD3M(IMOL,K),
     .                          ISCDEM(IMOL,K),IESTMM(IMOL,K),
     .                          IBGKM(IMOL,K)
            ELSE IF (ISRSM(IMOL).EQ.4) THEN
              READ (IUNIN,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                          ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                          ISCD3M(IMOL,K),ISCD4M(IMOL,K),
     .                          ISCDEM(IMOL,K),IESTMM(IMOL,K),
     .                          IBGKM(IMOL,K)
            ENDIF
            READ (IUNIN,6664) EELECM(IMOL,K),EBULKM(IMOL,K),
     .                        ESCD1M(IMOL,K),ESCD2M(IMOL,K),
     .                        FREACM(IMOL,K)
 4322     CONTINUE
          J=0
          DO K=1,NRCM(IMOL)
            IF (IREAC1(IREACM(IMOL,K)).GE.1) THEN
              J = J + 1
              IREACM(IMOL,J) = IREAC1(IREACM(IMOL,K))
              IBULKM(IMOL,J) = IBULKM(IMOL,K)
              ISCD1M(IMOL,J) = ISCD1M(IMOL,K)
              ISCD2M(IMOL,J) = ISCD2M(IMOL,K)
              ISCD3M(IMOL,J) = ISCD3M(IMOL,K)
              ISCD4M(IMOL,J) = ISCD4M(IMOL,K)
              ISCDEM(IMOL,J) = ISCDEM(IMOL,K)
              IESTMM(IMOL,J) = IESTMM(IMOL,K)
              IBGKM(IMOL,J) = IBGKM(IMOL,K)
              IF (B2IDEZ(ISCDEM(IMOL,J),5,5).eq.3) THEN
                IR = nint(EELECM(IMOL,K))
                EELECM(IMOL,J) = dble(IREAC1(IR))
              ELSE
                EELECM(IMOL,J) = EELECM(IMOL,K)
              END IF
              EBULKM(IMOL,J) = EBULKM(IMOL,K)
              ESCD1M(IMOL,J) = ESCD1M(IMOL,K)
              ESCD2M(IMOL,J) = ESCD2M(IMOL,K)
              FREACM(IMOL,J) = FREACM(IMOL,K)
            END IF
          END DO
          NRCM(IMOL) = J
4312   CONTINUE
*4C.   TEST IONS SPECIES CARDS, NIONI SPECIES
        READ (IUNIN,'(A72)') TXTION
  440   READ (IUNIN,*) NIONI
        allocate(TEXTI(NIONI))
        allocate(NMASSI(NIONI),NCHARI(NIONI),
     .           NPRTI(NIONI),NCHRGI(NIONI),
     .           ISRFI(NIONI),ISRTI(NIONI),ISRSI(NIONI),
     .           NRCI(NIONI),NHSTI(NIONI))
        allocate(NFOLI(NIONI),NGENI(NIONI))
        DO 441 IION=1,NIONI
          IF (streql(eirene_format,'juelich') .or.
     .        streql(eirene_format,'iter')) THEN
            READ (IUNIN,66666) I,TEXTI(IION),NMASSI(IION),NCHARI(IION),
     .                         NPRTI(IION),NCHRGI(IION),
     .                         ISRFI(IION),ISRTI(IION),ISRSI(IION),
     .                         NRCI(IION),NFOLI(IION),NGENI(IION),
     .                         NHSTI(IION),LKINDI(IION)
          ELSE IF (streql(eirene_format,'facelift') .or.
     .             streql(eirene_format,'new')) THEN
            READ (IUNIN,66666) I,TEXTI(IION),NMASSI(IION),NCHARI(IION),
     .                         NPRTI(IION),NCHRGI(IION),
     .                         ISRFI(IION),ISRTI(IION),ISRSI(IION),
     .                         NRCI(IION),NFOLI(IION),NGENI(IION)
            NHSTI(IION) = 0
            LKINDI(IION) = LIONSCL(IION)
          ELSE IF (streql(eirene_format,'old')) THEN
            READ (IUNIN,66666) I,TEXTI(IION),NMASSI(IION),NCHARI(IION),
     .                         NPRTI(IION),NCHRGI(IION),
     .                         ISRFI(IION),ISRTI(IION),ISRSI(IION),
     .                         NRCI(IION),NFOLI(IION),LKINDI(IION)
            NGENI(IION) = 0
            NHSTI(IION) = 0
          ENDIF
          DO 442 K=1,NRCI(IION)
            IF (ISRSI(IION).LT.3) THEN
              READ (IUNIN,6666) idum(0:6) ! IREACI(IION,K),IBULKI(IION,K),
                                          ! ISCD1I(IION,K),ISCD2I(IION,K),
                                          ! ISCDEI(IION,K),IESTMI(IION,K),
                                          ! IBGKI(IION,K)
            ELSE IF (ISRSI(IION).EQ.3) THEN
              READ (IUNIN,6666) idum(0:7)
            ELSE IF (ISRSI(IION).EQ.4) THEN
              READ (IUNIN,6666) idum(0:8)
            ENDIF
            READ (IUNIN,6664) rdum(0:4) ! EELECI(IION,K),EBULKI(IION,K),
                                        ! ESCD1I(IION,K),ESCD2I(IION,K),
                                        ! FREACI(IION,K)
442       CONTINUE
441     CONTINUE
        DO 4411 IION=1,NIONI
          backspace(IUNIN)
          DO 4421 K=1,NRCI(IION)
            backspace(IUNIN)
            backspace(IUNIN)
 4421     CONTINUE
 4411   CONTINUE
        NNREAC=max(1,maxval(NRCI(1:NIONI)))
        allocate(IREACI(NIONI,NNREAC),IBULKI(NIONI,NNREAC),
     .           ISCD1I(NIONI,NNREAC),ISCD2I(NIONI,NNREAC),
     .           ISCD3I(NIONI,NNREAC),ISCD4I(NIONI,NNREAC),
     .           ISCDEI(NIONI,NNREAC))
        allocate(IESTMI(NIONI,NNREAC),IBGKI(NIONI,NNREAC))
        allocate(EELECI(NIONI,NNREAC),EBULKI(NIONI,NNREAC),
     .           ESCD1I(NIONI,NNREAC),ESCD2I(NIONI,NNREAC),
     .           FREACI(NIONI,NNREAC))
        DO 4412 IION=1,NIONI
          READ (IUNIN,'(A72)') ZEILE
          DO 4422 K=1,NRCI(IION)
            IF (ISRSI(IION).LT.3) THEN
              READ (IUNIN,6666) IREACI(IION,K),IBULKI(IION,K),
     .                          ISCD1I(IION,K),ISCD2I(IION,K),
     .                          ISCDEI(IION,K),IESTMI(IION,K),
     .                          IBGKI(IION,K)
            ELSE IF (ISRSI(IION).EQ.3) THEN
              READ (IUNIN,6666) IREACI(IION,K),IBULKI(IION,K),
     .                          ISCD1I(IION,K),ISCD2I(IION,K),
     .                          ISCD3I(IION,K),
     .                          ISCDEI(IION,K),IESTMI(IION,K),
     .                          IBGKI(IION,K)
            ELSE IF (ISRSI(IION).EQ.4) THEN
              READ (IUNIN,6666) IREACI(IION,K),IBULKI(IION,K),
     .                          ISCD1I(IION,K),ISCD2I(IION,K),
     .                          ISCD3I(IION,K),ISCD4I(IION,K),
     .                          ISCDEI(IION,K),IESTMI(IION,K),
     .                          IBGKI(IION,K)
            ENDIF
            READ (IUNIN,6664) EELECI(IION,K),EBULKI(IION,K),
     .                        ESCD1I(IION,K),ESCD2I(IION,K),
     .                        FREACI(IION,K)
 4422     CONTINUE
          J=0
          DO K=1,NRCI(IION)
            IF (IREAC1(IREACI(IION,K)).GE.1) THEN
              J = J + 1
              IREACI(IION,J) = IREAC1(IREACI(IION,K))
              IBULKI(IION,J) = IBULKI(IION,K)
              ISCD1I(IION,J) = ISCD1I(IION,K)
              ISCD2I(IION,J) = ISCD2I(IION,K)
              ISCD3I(IION,J) = ISCD3I(IION,K)
              ISCD4I(IION,J) = ISCD4I(IION,K)
              ISCDEI(IION,J) = ISCDEI(IION,K)
              IESTMI(IION,J) = IESTMI(IION,K)
              IBGKI(IION,J) = IBGKI(IION,K)
              IF (B2IDEZ(ISCDEI(IION,J),5,5).eq.3) THEN
                IR = nint(EELECI(IION,K))
                EELECI(IION,J) = dble(IREAC1(IR))
              ELSE
                EELECI(IION,J) = EELECI(IION,K)
              END IF
              EBULKI(IION,J) = EBULKI(IION,K)
              ESCD1I(IION,J) = ESCD1I(IION,K)
              ESCD2I(IION,J) = ESCD2I(IION,K)
              FREACI(IION,J) = FREACI(IION,K)
            END IF
          END DO
          NRCI(IION) = J
 4412   CONTINUE
*4D.   NEUTRAL PHOTONS SPECIES CARDS, NPHOTI SPECIES
        READ (IUNIN,'(A80)') ZEILE
        IF (ZEILE(1:3).EQ.'***') THEN
          IREAD=1
          NPHOTI=0
          NCRM=0
          READ (ZEILE,'(A72)') TITLE5
          TXTPHOT=''
          GOTO 500
        ELSE
          READ (ZEILE,'(A72)') TXTPHOT
        END IF
        IREAD=0
  450   READ (IUNIN,*) NPHOTI
        allocate(TEXTPH(NPHOTI))
        allocate(ISRFPH(NPHOTI),ISRTPH(NPHOTI),ISRSPH(NPHOTI),
     .           NRCPH(NPHOTI))
        allocate(NFOLPH(NPHOTI),NGENPH(NPHOTI),NHSTPH(NPHOTI))
        DO 451 IPHOT=1,NPHOTI
          READ (IUNIN,66666) I,TEXTPH(IPHOT),idum(0:3),
     .                       ISRFPH(IPHOT),ISRTPH(IPHOT),ISRSPH(IPHOT),
     .                       NRCPH(IPHOT),NFOLPH(IPHOT),NGENPH(IPHOT),
     .                       NHSTPH(IPHOT)
          DO 452 K=1,NRCPH(IPHOT)
            IF (ISRSPH(IPHOT).LT.3) THEN
              READ (IUNIN,6666) idum(0:6) ! IREACPH(IPHOT,K),IBULKPH(IPHOT,K),
                                          ! ISCD1PH(IPHOT,K),ISCD2PH(IPHOT,K),
                                          ! ISCDEPH(IPHOT,K),IESTMPH(IPHOT,K),
                                          ! IBGKPH(IPHOT,K)
            ELSE IF (ISRSPH(IPHOT).EQ.3) THEN
              READ (IUNIN,6666) idum(0:7)
            ELSE IF (ISRSPH(IPHOT).EQ.4) THEN
              READ (IUNIN,6666) idum(0:8)
            ENDIF
            READ (IUNIN,6664) rdum(0:5) ! EELECPH(IPHOT,K),EBULKPH(IPHOT,K),
                                        ! ESCD1PH(IPHOT,K),ESCD2PH(IPHOT,K),
                                        ! FREACPH(IPHOT,K),FLDLMPH(IPHOT,K)
452       CONTINUE
451     CONTINUE
        DO 4511 IPHOT=1,NPHOTI
          backspace(IUNIN)
          DO 4521 K=1,NRCPH(IPHOT)
            backspace(IUNIN)
            backspace(IUNIN)
 4521     CONTINUE
 4511   CONTINUE
        NNREAC=max(1,maxval(NRCPH(1:NPHOTI)))
        allocate(IREACPH(NPHOTI,NNREAC),IBULKPH(NPHOTI,NNREAC),
     .           ISCD1PH(NPHOTI,NNREAC),ISCD2PH(NPHOTI,NNREAC),
     .           ISCD3PH(NPHOTI,NNREAC),ISCD4PH(NPHOTI,NNREAC),
     .           ISCDEPH(NPHOTI,NNREAC))
        allocate(IESTMPH(NPHOTI,NNREAC),IBGKPH(NPHOTI,NNREAC))
        allocate(EELECPH(NPHOTI,NNREAC),EBULKPH(NPHOTI,NNREAC),
     .           ESCD1PH(NPHOTI,NNREAC),ESCD2PH(NPHOTI,NNREAC),
     .           FREACPH(NPHOTI,NNREAC),FLDLMPH(NPHOTI,NNREAC))
        DO 4512 IPHOT=1,NPHOTI
          READ (IUNIN,'(A72)') ZEILE
          DO 4522 K=1,NRCPH(IPHOT)
            IF (ISRSPH(IPHOT).LT.3) THEN
              READ (IUNIN,6666) IREACPH(IPHOT,K),IBULKPH(IPHOT,K),
     .                          ISCD1PH(IPHOT,K),ISCD2PH(IPHOT,K),
     .                          ISCDEPH(IPHOT,K),IESTMPH(IPHOT,K),
     .                          IBGKPH(IPHOT,K)
            ELSE IF (ISRSPH(IPHOT).EQ.3) THEN
              READ (IUNIN,6666) IREACPH(IPHOT,K),IBULKPH(IPHOT,K),
     .                          ISCD1PH(IPHOT,K),ISCD2PH(IPHOT,K),
     .                          ISCD3PH(IPHOT,K),
     .                          ISCDEPH(IPHOT,K),IESTMPH(IPHOT,K),
     .                          IBGKPH(IPHOT,K)
            ELSE IF (ISRSPH(IPHOT).EQ.4) THEN
              READ (IUNIN,6666) IREACPH(IPHOT,K),IBULKPH(IPHOT,K),
     .                          ISCD1PH(IPHOT,K),ISCD2PH(IPHOT,K),
     .                          ISCD3PH(IPHOT,K),ISCD4PH(IPHOT,K),
     .                          ISCDEPH(IPHOT,K),IESTMPH(IPHOT,K),
     .                          IBGKPH(IPHOT,K)
            ENDIF
            READ (IUNIN,6664) EELECPH(IPHOT,K),EBULKPH(IPHOT,K),
     .                        ESCD1PH(IPHOT,K),ESCD2PH(IPHOT,K),
     .                        FREACPH(IPHOT,K),FLDLMPH(IPHOT,K)
 4522     CONTINUE
          J=0
          DO K=1,NRCPH(IPHOT)
            IF (IREAC1(IREACPH(IPHOT,K)).GE.1) THEN
              J = J + 1
              IREACPH(IPHOT,J) = IREAC1(IREACPH(IPHOT,K))
              IBULKPH(IPHOT,J) = IBULKPH(IPHOT,K)
              ISCD1PH(IPHOT,J) = ISCD1PH(IPHOT,K)
              ISCD2PH(IPHOT,J) = ISCD2PH(IPHOT,K)
              ISCD3PH(IPHOT,J) = ISCD3PH(IPHOT,K)
              ISCD4PH(IPHOT,J) = ISCD4PH(IPHOT,K)
              ISCDEPH(IPHOT,J) = ISCDEPH(IPHOT,K)
              IESTMPH(IPHOT,J) = IESTMPH(IPHOT,K)
              IBGKPH(IPHOT,J) = IBGKPH(IPHOT,K)
              IF (B2IDEZ(ISCDEPH(IPHOT,J),5,5).eq.3) THEN
                IR = nint(EELECPH(IPHOT,K))
                EELECPH(IPHOT,J) = dble(IREAC1(IR))
              ELSE
                EELECPH(IPHOT,J) = EELECPH(IPHOT,K)
              END IF
              EBULKPH(IPHOT,J) = EBULKPH(IPHOT,K)
              ESCD1PH(IPHOT,J) = ESCD1PH(IPHOT,K)
              ESCD2PH(IPHOT,J) = ESCD2PH(IPHOT,K)
              FREACPH(IPHOT,J) = FREACPH(IPHOT,K)
              FLDLMPH(IPHOT,J) = FLDLMPH(IPHOT,K)
            END IF
          END DO
          NRCPH(IPHOT) = J
 4512   CONTINUE
*4E. DATA FOR COLLISIONAL-RADIATIVE MODELS
        READ (IUNIN,'(A72)') ZEILE
        IF (ZEILE(1:3).EQ.'***') THEN
          READ (ZEILE,'(A72)') TITLE5
          IREAD=1
          GOTO 500
        ELSE
          READ (ZEILE,'(A72)') TXTCRM
        ENDIF
        READ (IUNIN,'(I6)') NCRM
        IREAD=0
        IF (NCRM.LE.0) GOTO 500
        allocate(NSCRM(NCRM),NIONREAC(NCRM),NRCMBREAC(NCRM),
     .           TYPCRM(NCRM))
        DO ICRM=1,NCRM
          READ(IUNIN,'(5I6)') I,NSCRM(ICRM),NIONREAC(ICRM),
     .                        NRCMBREAC(ICRM),TYPCRM(ICRM)
          DO K=1,NSCRM(ICRM)
            READ(IUNIN,'(4I6)') idum(0:3)
          END DO
        END DO
        NNREAC=MAX(1,maxval(NSCRM(1:NCRM)))
        DO ICRM=1,NCRM
          backspace(IUNIN)
          DO K=1,NSCRM(ICRM)
            backspace(IUNIN)
          END DO
        END DO
        allocate(SCRM(NNREAC,NCRM),NNCRM(NNREAC,NCRM),
     .          PHCRM(NNREAC,NCRM))
        DO ICRM=1,NCRM
          READ(IUNIN,'(5I6)') idum(0:4)
          DO K=1,NSCRM(ICRM)
            READ(IUNIN,'(4I6)') I,SCRM(K,ICRM),
     r                          NNCRM(K,ICRM),PHCRM(K,ICRM)
          END DO
        END DO
*** 5. DATA FOR PLASMA BACKGROUND
  500   CONTINUE
        IF (IREAD.EQ.0) READ (IUNIN,'(A72)') TITLE5
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
        allocate(NHSTP(NPLSI))
        allocate(DENSLIM(NPLSI))
        allocate(CDENMODEL(NPLSI))
        allocate(NRE(NPLSI))
        allocate(TDMBOLTZ(NPLSI))
        allocate(TDMDELTAE(NPLSI))
        allocate(TDM_CORONA(NPLSI))
        DO 511 IPLS=1,NPLSI
          READ (IUNIN,66666) I,TEXTP(IPLS),NMASSP(IPLS),NCHARP(IPLS),
     .                                     NPRTP(IPLS),NCHRGP(IPLS),
     .                                     ISRFP(IPLS),ISRTP(IPLS),
     .                                     ISRSP(IPLS),NRCP(IPLS),
     .                                     DENSLIM(IPLS),idum(0),
     .                                     NHSTP(IPLS),idum(2),
     .                                     CDENMODEL(IPLS),NRE(IPLS)
          CALL EIRENE_UPPERCASE (CDENMODEL(IPLS))
          DO 512 K=1,NRCP(IPLS)
            IF (ISRSP(IPLS).LT.3) THEN
              READ (IUNIN,6666) idum(0:4)
            ELSE IF (ISRSP(IPLS).EQ.3) THEN
              READ (IUNIN,6666) idum(0:5)
            ELSE IF (ISRSP(IPLS).EQ.4) THEN
              READ (IUNIN,6666) idum(0:6)
            ENDIF
            READ (IUNIN,6664) rdum(0:4)
512       CONTINUE
          IF (CDENMODEL(IPLS).EQ.'FORT.13   ' .OR.
     &        CDENMODEL(IPLS).EQ.'BOLTZMANN ' .OR.
     &        CDENMODEL(IPLS).EQ.'CORONA    ') THEN
            READ (IUNIN,'(A80)') ZEILE
          ELSE IF (CDENMODEL(IPLS).EQ.'COLRAD    ') THEN
            DO K=1,NRE(IPLS)
              READ (IUNIN,'(A80)') ZEILE
            ENDDO
          ENDIF
511     CONTINUE
        DO 5111 IPLS=1,NPLSI
          backspace(IUNIN)
          DO 5121 K=1,NRCP(IPLS)
            backspace(IUNIN)
            backspace(IUNIN)
 5121     CONTINUE
          IF (CDENMODEL(IPLS).EQ.'FORT.13   ' .OR.
     &        CDENMODEL(IPLS).EQ.'BOLTZMANN ' .OR.
     &        CDENMODEL(IPLS).EQ.'CORONA    ') THEN
            backspace(IUNIN)
          ELSE IF (CDENMODEL(IPLS).EQ.'COLRAD    ') THEN
            DO K=1,NRE(IPLS)
              backspace(IUNIN)
            ENDDO
          ENDIF
 5111   CONTINUE
        NNREAC=max(1,maxval(NRCP(1:NPLSI)))
        allocate(IREACP(NPLSI,NNREAC),IBULKP(NPLSI,NNREAC),
     .           ISCD1P(NPLSI,NNREAC),ISCD2P(NPLSI,NNREAC),
     .           ISCD3P(NPLSI,NNREAC),ISCD4P(NPLSI,NNREAC),
     .           ISCDEP(NPLSI,NNREAC))
        allocate(EELECP(NPLSI,NNREAC),EBULKP(NPLSI,NNREAC),
     .           ESCD1P(NPLSI,NNREAC),ESCD2P(NPLSI,NNREAC),
     .           FREACP(NPLSI,NNREAC))
        NNREAC=max(1,maxval(NRE(1:NPLSI)))
        allocate(TDMISP(NNREAC,NPLSI))
        allocate(TDMFNAME(NNREAC,NPLSI), TDMH2(NNREAC,NPLSI),
     .           TDMREACTION(NNREAC,NPLSI), TDMCR(NNREAC,NPLSI))
        DO 5112 IPLS=1,NPLSI
          READ (IUNIN,'(A72)') ZEILE
          DO 5122 K=1,NRCP(IPLS)
            IF (ISRSP(IPLS).LT.3) THEN
              READ (IUNIN,6666) IREACP(IPLS,K),IBULKP(IPLS,K),
     .                          ISCD1P(IPLS,K),ISCD2P(IPLS,K),
     .                          ISCDEP(IPLS,K)
            ELSE IF (ISRSP(IPLS).EQ.3) THEN
              READ (IUNIN,6666) IREACP(IPLS,K),IBULKP(IPLS,K),
     .                          ISCD1P(IPLS,K),ISCD2P(IPLS,K),
     .                          ISCD3P(IPLS,K),ISCDEP(IPLS,K)
            ELSE IF (ISRSP(IPLS).EQ.4) THEN
              READ (IUNIN,6666) IREACP(IPLS,K),IBULKP(IPLS,K),
     .                          ISCD1P(IPLS,K),ISCD2P(IPLS,K),
     .                          ISCD3P(IPLS,K),ISCD4P(IPLS,K),
     .                          ISCDEP(IPLS,K)
            ENDIF
            READ (IUNIN,6664) EELECP(IPLS,K),EBULKP(IPLS,K),
     .                        ESCD1P(IPLS,K),ESCD2P(IPLS,K),
     .                        FREACP(IPLS,K)
 5122     CONTINUE
          J=0
          DO K=1,NRCP(IPLS)
            IF (IREAC1(IREACP(IPLS,K)).GE.1) THEN
              J = J + 1
              IREACP(IPLS,J) = IREAC1(IREACP(IPLS,K))
              IBULKP(IPLS,J) = IBULKP(IPLS,K)
              ISCD1P(IPLS,J) = ISCD1P(IPLS,K)
              ISCD2P(IPLS,J) = ISCD2P(IPLS,K)
              ISCD3P(IPLS,J) = ISCD3P(IPLS,K)
              ISCD4P(IPLS,J) = ISCD4P(IPLS,K)
              ISCDEP(IPLS,J) = ISCDEP(IPLS,K)
              IF (B2IDEZ(ISCDEP(IPLS,J),5,5).eq.3) THEN
                IR = nint(EELECP(IPLS,K))
                EELECP(IPLS,J) = dble(IREAC1(IR))
              ELSE
                EELECP(IPLS,J) = EELECP(IPLS,K)
              END IF
              EBULKP(IPLS,J) = EBULKP(IPLS,K)
              ESCD1P(IPLS,J) = ESCD1P(IPLS,K)
              ESCD2P(IPLS,J) = ESCD2P(IPLS,K)
              FREACP(IPLS,J) = FREACP(IPLS,K)
            END IF
          END DO
          NRCP(IPLS) = J
          SELECT CASE (CDENMODEL(IPLS))
          CASE ('FORT.13   ')
            READ (IUNIN,6666) TDMISP(1,IPLS)
          CASE ('BOLTZMANN ')
            READ (IUNIN,'(I6,6X,2E12.4)') TDMISP(1,IPLS),
     .           TDMBOLTZ(IPLS), TDMDELTAE(IPLS)
          CASE ('CORONA    ')
            READ (IUNIN,'(I6,1X,A6,1X,A4,A9,A3,E12.4)')
     .           TDMISP(1,IPLS),TDMFNAME(1,IPLS),
     .           TDMH2(1,IPLS),TDMREACTION(1,IPLS),
     .           TDMCR(1,IPLS),TDM_CORONA(IPLS)
          CASE ('COLRAD    ')
            DO I=1, NRE(IPLS)
              READ (IUNIN,'(I6,1X,A6,1X,A4,A9,A3,E12.4)')
     .             TDMISP(I,IPLS),TDMFNAME(I,IPLS),
     .             TDMH2(I,IPLS),TDMREACTION(I,IPLS),
     .             TDMCR(I,IPLS)
            END DO
          CASE DEFAULT
            CONTINUE
          END SELECT
 5112   CONTINUE
cxpb Now we have read the species data and can identify the bundles
        DO 5113 IPLS=1,NPLSI
          found = ipls.gt.(ns-natmi)
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
          end if
 5113   CONTINUE
        NPLSI1 = (ns1-natmi1)+(NPLSI-(ns-natmi))
        allocate(NCHARP1(NPLSI1))
        allocate(NMASSP1(NPLSI1))
        allocate(NCHRGP1(NPLSI1))
        allocate(ISP(NPLSI1))
        allocate(e_bundle_start(nplsi1))
        allocate(e_bundle_end(nplsi1))
        allocate(e_bundle_contr(nplsi1,npls))
        ib = 0
        do IATM = 1, NATMI1
          if (IATM.eq.NATMI1) then
            iend = ns1
          else
            iend = eb2spcr1(IATM+1)
          endif
          do iss = eb2spcr1(IATM)+1, iend-1
            ib = ib + 1
            NMASSP1(ib) = nint(am1(iss))
            NCHARP1(ib) = nint(zn1(iss))
            NCHRGP1(ib) = nint((zamin1(iss)+zamax1(iss))/2.0_R8)
          enddo
        end do
        e_bundle_start = 0
        e_bundle_end = 0
        e_bundle_contr = 0.0_R8
        ipls1 = 0
        ISP = 0
        ipls = 0
        iplss = 0
        do IATM = 1, NATMI
          if (iatm1(iatm).eq.0) then
            do while (NMASSP(ipls+1).eq.NMASSA(iatm) .and.
     .                NCHARP(ipls+1).eq.NCHARA(iatm))
               ipls = ipls + 1
            end do
            cycle
          endif
          is = eb2spcr(iatm)+1
          if (iatm.eq.natmi) then
            js = ns-1
          else
            js = eb2spcr(iatm+1)-1
          endif
          iss = eb2spcr1(iatm1(iatm))+1
          jss = iss
          if (jss.lt.ns1-1) then
            do while (jss.le.ns1-2 .and.
     .       am1(min(jss+1,ns1-1)).eq.am1(iss) .and.
     .       zn1(min(jss+1,ns1-1)).eq.zn1(iss))
              jss = jss + 1
            end do
          end if
          do i = is, js
            ipls = ipls+1
            do j = iss, jss
              if (iplss.lt.nplsi1) then
                if (NMASSP1(iplss+1).eq.NMASSP(ipls) .and.
     .              NCHARP1(iplss+1).eq.NCHARP(ipls)) then
                  zmin1 = zamin(i)
                  zmax1 = zamax(i)
                  zmin2 = zamin1(j)
                  zmax2 = zamax1(j)
                  charges_match =
     .             (zmin1.ge.zmin2 .and. zmin1.le.zmax2) .or.
     .             (zmax1.ge.zmin2 .and. zmax1.le.zmax2) .or.
     .             (zmin1.le.zmin2 .and. zmax1.ge.zmax2)
                  if (charges_match) then
cxpb  We must only increase the index number of the new plasma
cxpb  species if the charge bundle from the old plasma species
cxpb  is not fully included in the new species and is not the
cxpb  first one
                    if (.not.(zmin1.gt.zmin2 .and. zmax1.le.zmax2))
     >               iplss = iplss + 1
                    ipls1(ipls) = iplss
                    if (isp(iplss).eq.0) isp(iplss) = ipls
                    if (e_bundle_start(iplss).eq.0) then
                      e_bundle_start(iplss) = ipls
                      e_bundle_end(iplss) = ipls
                    else
                      e_bundle_end(iplss) = ipls
                    endif
                    e_bundle_contr(iplss,ipls) =
     .               (min(zmax1,zmax2)-max(zmin1,zmin2)+1.0_R8)/
     .               (zmax1-zmin1+1.0_R8)
                  endif
                endif
              else
                if (NMASSP1(iplss).eq.NMASSP(ipls).and.
     .              NCHARP1(iplss).eq.NCHARP(ipls)) then
                  zmin1 = zamin(i)
                  zmax1 = zamax(i)
                  zmin2 = zamin1(j)
                  zmax2 = zamax1(j)
                  charges_match =
     .             (zmin1.ge.zmin2 .and. zmin1.le.zmax2) .or.
     .             (zmax1.ge.zmin2 .and. zmax1.le.zmax2) .or.
     .             (zmin1.le.zmin2 .and. zmax1.ge.zmax2)
                  if (charges_match) then
                    ipls1(ipls) = iplss
                    if (isp(iplss).eq.0) isp(iplss) = ipls
                    if (e_bundle_start(iplss).eq.0) then
                      e_bundle_start(iplss) = ipls
                      e_bundle_end(iplss) = ipls
                    else
                      e_bundle_end(iplss) = ipls
                    endif
                    e_bundle_contr(iplss,ipls) =
     .               (min(zmax1,zmax2)-max(zmin1,zmin2)+1.0_R8)/
     .               (zmax1-zmin1+1.0_R8)
                  endif
                endif
              endif
            enddo
          enddo
        enddo
        do IPLS = ns-natmi+1, NPLSI
          ib = (ns1-natmi1) + IPLS-(ns-natmi)
          NMASSP1(ib) = NMASSP(IPLS)
          NCHARP1(ib) = NCHARP(IPLS)
          NCHRGP1(ib) = NCHRGP(IPLS)
          ipls1(ipls) = ib
          isp(ib) = ipls
          e_bundle_start(ib) = ipls
          e_bundle_end(ib) = ipls
          e_bundle_contr(ib,ipls) = 1.0_R8
        end do
        NPLSM = max(NPLSI,NPLSI1)
cxpb Must now convert the bulk ion data according to the new bundling
        allocate(TEXTP1(NPLSI1))
        allocate(NPRTP1(NPLSI1))
        allocate(ISRFP1(NPLSI1))
        allocate(ISRTP1(NPLSI1))
        allocate(ISRSP1(NPLSI1))
        allocate(NRCP1(NPLSI1))
        allocate(NHSTP1(NPLSI1))
        allocate(DENSLIM1(NPLSI1))
        allocate(CDENMODEL1(NPLSI1))
        allocate(NRE1(NPLSI1))
        allocate(TDMBOLTZ1(NPLSI1))
        allocate(TDMDELTAE1(NPLSI1))
        allocate(TDM_CORONA1(NPLSI1))
        allocate(done(NPLSI))
        done = .false.
        DO IB = 1, NPLSI1
          IS = ISP(IB)
          IF (IS.NE.0) THEN
            TEXTP1(IB) = ' '
            IF (e_bundle_start(IB).eq.e_bundle_end(IB)) then
              if (.not.done(IS)) then
                TEXTP1(IB) = TEXTP(e_bundle_start(IB))
              else
                len = lnblnk(TEXTP(e_bundle_start(IB)))
                if (NCHRGP(e_bundle_start(IB)).lt.10) then
                  write(TEXTP1(IB)(1:len-2),'(A)')
     .             TEXTP(e_bundle_start(IB))(1:len-2)
                else
                  write(TEXTP1(IB)(1:len-3),'(A)')
     .             TEXTP(e_bundle_start(IB))(1:len-3)
                endif
                if (NCHRGP1(IB).lt.10) then
                  write(TEXTP1(IB)(len-1:len-1),'(i1)') NCHRGP1(IB)
                else
                  write(TEXTP1(IB)(len-2:len-1),'(i2)') NCHRGP1(IB)
                endif
                write(TEXTP1(IB)(len:len),'(a1)') '+'
              endif
            ELSE
              TEXTP1(IB) =
     . TEXTP(e_bundle_start(IB))(1:lnblnk(TEXTP(e_bundle_start(IB))))
     . //'_'//
     . TEXTP(e_bundle_end(IB))(1:lnblnk(TEXTP(e_bundle_end(IB))))
            END IF
            NPRTP1(IB) = 0
            NHSTP1(IB) = -1
            NRCP1(IB) = 0
            do ipls = e_bundle_start(IB), e_bundle_end(IB)
              NPRTP1(IB) = NPRTP1(IB) + NPRTP(ipls)
              NRCP1(IB) = NRCP1(IB) + NRCP(ipls)
              IF (NHSTP(ipls).EQ.0) NHSTP1(IB) = 0
            end do
            NPRTP1(IB) =
     .       NPRTP1(IB) / (e_bundle_end(IB)-e_bundle_start(IB)+1)
            NRE1(IB) = 0
            IF (ISRFP(e_bundle_end(IB)).GT.0.AND.
     .          ISRFP(e_bundle_end(IB)).LE.NATMI) THEN
              ISRFP1(IB) = IATM1(ISRFP(e_bundle_end(IB)))
            ELSE
              ISRFP1(IB) = 0
            ENDIF
            IF (ABS(ISRTP(e_bundle_end(IB))).GT.0.AND.
     .          ABS(ISRTP(e_bundle_end(IB))).LE.NATMI) THEN
              ISRTP1(IB) =
     .         SIGN(IATM1(ABS(ISRTP(e_bundle_end(IB)))),
     .                        ISRTP(e_bundle_end(IB)))
            ELSE
              ISRTP1(IB) = 0
            ENDIF
            ISRSP1(IB) = ISRSP(e_bundle_end(IB))
            DENSLIM1(IB) = DENSLIM(e_bundle_end(IB))
            CDENMODEL1(IB) = CDENMODEL(e_bundle_end(IB))
            DO IPLS = e_bundle_start(IB), e_bundle_end(IB)
              IF (CDENMODEL1(IB)(1:6).EQ.'COLRAD') THEN
                NRE1(IB) = NRE1(IB) + NRE(IPLS)
              ELSE
                NRE1(IB) = 1
              END IF
            END DO
            done(IS) = .true.
          END IF
        END DO
        NREAC1 = max(1,maxval(NRE1(1:NPLSI1)))
        allocate(TDMISP1(NREAC1,NPLSI1))
        allocate(TDMFNAME1(NREAC1,NPLSI1))
        allocate(TDMH21(NREAC1,NPLSI1))
        allocate(TDMREACTION1(NREAC1,NPLSI1))
        allocate(TDMCR1(NREAC1,NPLSI1))
        DO IB = 1, NPLSI1
          IS = ISP(IB)
          IF (IS.NE.0) THEN
            SELECT CASE (CDENMODEL1(IB))
            CASE ('FORT.13   ')
              TDMISP1(1,IB) = ipls1(TDMISP(1,e_bundle_end(IB)))
            CASE ('BOLTZMANN ')
              TDMISP1(1,IB) = ipls1(TDMISP(1,e_bundle_end(IB)))
              TDMBOLTZ1(IB) = TDMBOLTZ(e_bundle_end(IB))
              TDMDELTAE1(IB) = TDMDELTAE(e_bundle_end(IB))
            CASE ('CORONA    ')
              TDMISP1(1,IB) = ipls1(TDMISP(1,e_bundle_end(IB)))
              TDMFNAME1(1,IB) = TDMFNAME(1,e_bundle_end(IB))
              TDMH21(1,IB) = TDMH2(1,e_bundle_end(IB))
              TDMREACTION1(1,IB) = TDMREACTION(1,e_bundle_end(IB))
              TDMCR1(1,IB) = TDMCR(1,e_bundle_end(IB))
              TDM_CORONA1(IB) = TDM_CORONA(e_bundle_end(IB))
            CASE ('COLRAD    ')
              I1 = 0
              do ipls = e_bundle_start(IB), e_bundle_end(IB)
                DO I = 1, NRE(ipls)
                  I1 = I1 + 1
                  TDMISP1(I1,IB) = ipls1(TDMISP(I,IPLS))
                  TDMFNAME1(I1,IB) = TDMFNAME(I,IPLS)
                  TDMH21(I1,IB) = TDMH2(I,IPLS)
                  TDMREACTION1(I1,IB) = TDMREACTION(I,IPLS)
                  TDMCR1(I1,IB) = TDMCR(I,IPLS)
                END DO
              end do
            CASE DEFAULT
              CONTINUE
            END SELECT
          END IF
        END DO
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
            IF (ISRSA(IATM).GE.3) THEN
              M = B2IDEZ(ISCD3A(IATM,K),1,2)
              IF (M.EQ.1) THEN   ! we have an atom
                L = B2IDEZ(ISCD3A(IATM,K),2,3)
                ICRM = B2IDEZ(ISCD3A(IATM,K),3,3)
                ISCD3A(IATM,K) = M + 10*L + 100*iatm1(icrm)
              ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
                L = B2IDEZ(ISCD3A(IATM,K),2,3)
                IPLS = B2IDEZ(ISCD3A(IATM,K),3,3)
                ISCD3A(IATM,K) = M + 10*L + 100*ipls1(ipls)
              END IF
            END IF
            IF (ISRSA(IATM).GE.4) THEN
              M = B2IDEZ(ISCD4A(IATM,K),1,2)
              IF (M.EQ.1) THEN   ! we have an atom
                L = B2IDEZ(ISCD4A(IATM,K),2,3)
                ICRM = B2IDEZ(ISCD4A(IATM,K),3,3)
                ISCD4A(IATM,K) = M + 10*L + 100*iatm1(icrm)
              ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
                L = B2IDEZ(ISCD4A(IATM,K),2,3)
                IPLS = B2IDEZ(ISCD4A(IATM,K),3,3)
                ISCD4A(IATM,K) = M + 10*L + 100*ipls1(ipls)
              END IF
            END IF
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
            IF (ISRSM(IMOL).GE.3) THEN
              M = B2IDEZ(ISCD3M(IMOL,K),1,2)
              IF (M.EQ.1) THEN   ! we have an atom
                L = B2IDEZ(ISCD3M(IMOL,K),2,3)
                IATM = B2IDEZ(ISCD3M(IMOL,K),3,3)
                ISCD3M(IMOL,K) = M + 10*L + 100*iatm1(iatm)
              ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
                L = B2IDEZ(ISCD3M(IMOL,K),2,3)
                IPLS = B2IDEZ(ISCD3M(IMOL,K),3,3)
                ISCD3M(IMOL,K) = M + 10*L + 100*ipls1(ipls)
              END IF
            END IF
            IF (ISRSM(IMOL).GE.4) THEN
              M = B2IDEZ(ISCD4M(IMOL,K),1,2)
              IF (M.EQ.1) THEN   ! we have an atom
                L = B2IDEZ(ISCD4M(IMOL,K),2,3)
                IATM = B2IDEZ(ISCD4M(IMOL,K),3,3)
                ISCD4M(IMOL,K) = M + 10*L + 100*iatm1(iatm)
              ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
                L = B2IDEZ(ISCD4M(IMOL,K),2,3)
                IPLS = B2IDEZ(ISCD4M(IMOL,K),3,3)
                ISCD4M(IMOL,K) = M + 10*L + 100*ipls1(ipls)
              END IF
            END IF
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
            IF (ISRSI(IION).GE.3) THEN
              M = B2IDEZ(ISCD3I(IION,K),1,2)
              IF (M.EQ.1) THEN   ! we have an atom
                L = B2IDEZ(ISCD3I(IION,K),2,3)
                IATM = B2IDEZ(ISCD3I(IION,K),3,3)
                ISCD3I(IION,K) = M + 10*L + 100*iatm1(iatm)
              ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
                L = B2IDEZ(ISCD3I(IION,K),2,3)
                IPLS = B2IDEZ(ISCD3I(IION,K),3,3)
                ISCD3I(IION,K) = M + 10*L + 100*ipls1(ipls)
              END IF
            ENDIF
            IF (ISRSI(IION).GE.4) THEN
              M = B2IDEZ(ISCD4I(IION,K),1,2)
              IF (M.EQ.1) THEN   ! we have an atom
                L = B2IDEZ(ISCD4I(IION,K),2,3)
                IATM = B2IDEZ(ISCD4I(IION,K),3,3)
                ISCD4I(IION,K) = M + 10*L + 100*iatm1(iatm)
              ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
                L = B2IDEZ(ISCD4I(IION,K),2,3)
                IPLS = B2IDEZ(ISCD4I(IION,K),3,3)
                ISCD4I(IION,K) = M + 10*L + 100*ipls1(ipls)
              END IF
            ENDIF
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
          END DO
        END DO
        DO IPHOT = 1, NPHOTI
          DO K = 1, NRCPH(IPHOT)
            M = B2IDEZ(IBULKPH(IPHOT,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(IBULKPH(IPHOT,K),2,3)
              IATM = B2IDEZ(IBULKPH(IPHOT,K),3,3)
              IBULKPH(IPHOT,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(IBULKPH(IPHOT,K),2,3)
              IPLS = B2IDEZ(IBULKPH(IPHOT,K),3,3)
              IBULKPH(IPHOT,K) = M + 10*L + 100*ipls1(ipls)
            END IF
            M = B2IDEZ(ISCD1PH(IPHOT,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(ISCD1PH(IPHOT,K),2,3)
              IATM = B2IDEZ(ISCD1PH(IPHOT,K),3,3)
              ISCD1PH(IPHOT,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(ISCD1PH(IPHOT,K),2,3)
              IPLS = B2IDEZ(ISCD1PH(IPHOT,K),3,3)
              ISCD1PH(IPHOT,K) = M + 10*L + 100*ipls1(ipls)
            END IF
            M = B2IDEZ(ISCD2PH(IPHOT,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(ISCD2PH(IPHOT,K),2,3)
              IATM = B2IDEZ(ISCD2PH(IPHOT,K),3,3)
              ISCD2PH(IPHOT,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(ISCD2PH(IPHOT,K),2,3)
              IPLS = B2IDEZ(ISCD2PH(IPHOT,K),3,3)
              ISCD2PH(IPHOT,K) = M + 10*L + 100*ipls1(ipls)
            END IF
            IF (ISRSPH(IPHOT).GE.3) THEN
              M = B2IDEZ(ISCD3PH(IPHOT,K),1,2)
              IF (M.EQ.1) THEN   ! we have an atom
                L = B2IDEZ(ISCD3PH(IPHOT,K),2,3)
                IATM = B2IDEZ(ISCD3PH(IPHOT,K),3,3)
                ISCD3PH(IPHOT,K) = M + 10*L + 100*iatm1(iatm)
              ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
                L = B2IDEZ(ISCD3PH(IPHOT,K),2,3)
                IPLS = B2IDEZ(ISCD3PH(IPHOT,K),3,3)
                ISCD3PH(IPHOT,K) = M + 10*L + 100*ipls1(ipls)
              END IF
            ENDIF
            IF (ISRSPH(IPHOT).GE.4) THEN
              M = B2IDEZ(ISCD4PH(IPHOT,K),1,2)
              IF (M.EQ.1) THEN   ! we have an atom
                L = B2IDEZ(ISCD4PH(IPHOT,K),2,3)
                IATM = B2IDEZ(ISCD4PH(IPHOT,K),3,3)
                ISCD4PH(IPHOT,K) = M + 10*L + 100*iatm1(iatm)
              ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
                L = B2IDEZ(ISCD4PH(IPHOT,K),2,3)
                IPLS = B2IDEZ(ISCD4PH(IPHOT,K),3,3)
                ISCD4PH(IPHOT,K) = M + 10*L + 100*ipls1(ipls)
              END IF
            ENDIF
            M = B2IDEZ(IBGKPH(IPHOT,K),1,3)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(IBGKPH(IPHOT,K),2,3)
              IATM = B2IDEZ(IBGKPH(IPHOT,K),3,3)
              IBGKPH(IPHOT,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(IBGKPH(IPHOT,K),2,3)
              IPLS = B2IDEZ(IBGKPH(IPHOT,K),3,3)
              IBGKPH(IPHOT,K) = M + 10*L + 100*ipls1(ipls)
            END IF
          END DO
        END DO
        DO IPLS = 1, NPLSI
         IF (IPLS1(IPLS).ne.0) then
          DO K = 1, NRCP(IPLS)
            M = B2IDEZ(IBULKP(IPLS,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(IBULKP(IPLS,K),2,3)
              IATM = B2IDEZ(IBULKP(IPLS,K),3,3)
              IBULKP(IPLS,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(IBULKP(IPLS,K),2,3)
              ICRM = B2IDEZ(IBULKP(IPLS,K),3,3)
              IBULKP(IPLS,K) = M + 10*L + 100*ipls1(icrm)
            END IF
            M = B2IDEZ(ISCD1P(IPLS,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(ISCD1P(IPLS,K),2,3)
              IATM = B2IDEZ(ISCD1P(IPLS,K),3,3)
              ISCD1P(IPLS,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(ISCD1P(IPLS,K),2,3)
              ICRM = B2IDEZ(ISCD1P(IPLS,K),3,3)
              ISCD1P(IPLS,K) = M + 10*L + 100*ipls1(icrm)
            END IF
            M = B2IDEZ(ISCD2P(IPLS,K),1,2)
            IF (M.EQ.1) THEN   ! we have an atom
              L = B2IDEZ(ISCD2P(IPLS,K),2,3)
              IATM = B2IDEZ(ISCD2P(IPLS,K),3,3)
              ISCD2P(IPLS,K) = M + 10*L + 100*iatm1(iatm)
            ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
              L = B2IDEZ(ISCD2P(IPLS,K),2,3)
              ICRM = B2IDEZ(ISCD2P(IPLS,K),3,3)
              ISCD2P(IPLS,K) = M + 10*L + 100*ipls1(icrm)
            END IF
            IF (ISRSP(IPLS).GE.3) THEN
              M = B2IDEZ(ISCD3P(IPLS,K),1,2)
              IF (M.EQ.1) THEN   ! we have an atom
                L = B2IDEZ(ISCD3P(IPLS,K),2,3)
                IATM = B2IDEZ(ISCD3P(IPLS,K),3,3)
                ISCD3P(IPLS,K) = M + 10*L + 100*iatm1(iatm)
              ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
                L = B2IDEZ(ISCD3P(IPLS,K),2,3)
                ICRM = B2IDEZ(ISCD3P(IPLS,K),3,3)
                ISCD3P(IPLS,K) = M + 10*L + 100*ipls1(icrm)
              END IF
            END IF
            IF (ISRSP(IPLS).GE.4) THEN
              M = B2IDEZ(ISCD4P(IPLS,K),1,2)
              IF (M.EQ.1) THEN   ! we have an atom
                L = B2IDEZ(ISCD4P(IPLS,K),2,3)
                IATM = B2IDEZ(ISCD4P(IPLS,K),3,3)
                ISCD4P(IPLS,K) = M + 10*L + 100*iatm1(iatm)
              ELSE IF (M.EQ.4) THEN   ! we have a bulk ion
                L = B2IDEZ(ISCD4P(IPLS,K),2,3)
                ICRM = B2IDEZ(ISCD4P(IPLS,K),3,3)
                ISCD4P(IPLS,K) = M + 10*L + 100*ipls1(icrm)
              END IF
            END IF
          END DO
         END IF
        END DO
cxpb We now rewrite the modified block 4
*4A.   NEUTRAL ATOMS SPECIES CARDS, NATMI1 SPECIES
        WRITE(IUNOT,6666) NATMI1
        DO 4213 IATM=1,NATMI
         IF (IATM1(IATM).NE.0) THEN
          IF (ISRFA(IATM).GT.0 .AND.ISRFA(IATM).LE.NATMI) THEN
            ISRFA(IATM) = IATM1(ISRFA(IATM))
          ELSE
            ISRFA(IATM) = 0
          ENDIF
          IF (ABS(ISRTA(IATM)).GT.0.AND.ABS(ISRTA(IATM)).LE.NATMI) THEN
            ISRTA(IATM) = SIGN(IATM1(ABS(ISRTA(IATM))),ISRTA(IATM))
          ELSE
            ISRTA(IATM) = 0
          ENDIF
          WRITE(IUNOT,66666) IATM1(IATM),
     .                         TEXTA(IATM),NMASSA(IATM),
     .                         NCHARA(IATM),NPRTA(IATM),0,
     .                         ISRFA(IATM),ISRTA(IATM),ISRSA(IATM),
     .                         NRCA(IATM),NFOLA(IATM),NGENA(IATM),
     .                         NHSTA(IATM)
          DO 4223 K=1,NRCA(IATM)
            IF (ISRSA(IATM).LT.3) THEN
              WRITE(IUNOT,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                          ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                          ISCDEA(IATM,K),IESTMA(IATM,K),
     .                          IBGKA(IATM,K)
            ELSE IF (ISRSA(IATM).EQ.3) THEN
              WRITE(IUNOT,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                          ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                          ISCD3A(IATM,K),
     .                          ISCDEA(IATM,K),IESTMA(IATM,K),
     .                          IBGKA(IATM,K)
            ELSE IF (ISRSA(IATM).EQ.4) THEN
              WRITE(IUNOT,6666) IREACA(IATM,K),IBULKA(IATM,K),
     .                          ISCD1A(IATM,K),ISCD2A(IATM,K),
     .                          ISCD3A(IATM,K),ISCD4A(IATM,K),
     .                          ISCDEA(IATM,K),IESTMA(IATM,K),
     .                          IBGKA(IATM,K)
            ENDIF
            WRITE(IUNOT,6664) EELECA(IATM,K),EBULKA(IATM,K),
     .                        ESCD1A(IATM,K),ESCD2A(IATM,K),
     .                        FREACA(IATM,K),FLDLMA(IATM,K)
 4223     CONTINUE
         END IF
 4213   CONTINUE
*4B.   NEUTRAL MOLECULE SPECIES CARDS, NMOLI SPECIES
        WRITE(IUNOT,'(A72)') TXTMOL
        WRITE(IUNOT,6666) NMOLI
        DO 4313 IMOL=1,NMOLI
          IF (ISRFM(IMOL).GT.0 .AND. ISRFM(IMOL).LE.NATMI) THEN
            ISRFM(IMOL) = IATM1(ISRFM(IMOL))
          ELSE
            ISRFM(IMOL) = 0
          ENDIF
          IF (ABS(ISRTM(IMOL)).GT.0.AND.ABS(ISRTM(IMOL)).LE.NATMI) THEN
            ISRTM(IMOL) = SIGN(IATM1(ABS(ISRTM(IMOL))),ISRTM(IMOL))
          ELSE IF (ABS(ISRTM(IMOL)).GT.NATMI.AND.
     .             ABS(ISRTM(IMOL)).LE.NATMI+NMOLI) THEN
            ISRTM(IMOL) =
     .       SIGN(ABS(ISRTM(IMOL))+NATMI1-NATMI,ISRTM(IMOL))
          ELSE
            ISRTM(IMOL) = 0
          ENDIF
          IF (LKINDM(IMOL).GT.0 .AND. LKINDM(IMOL).LE.NATMI) THEN
            LKINDM(IMOL) = IATM1(LKINDM(IMOL))
          ELSE
            LKINDM(IMOL) = 0
          ENDIF
          WRITE(IUNOT,66666) IMOL,TEXTM(IMOL),NMASSM(IMOL),NCHARM(IMOL),
     .                         NPRTM(IMOL),0,
     .                         ISRFM(IMOL),ISRTM(IMOL),ISRSM(IMOL),
     .                         NRCM(IMOL),NFOLM(IMOL),NGENM(IMOL),
     .                         NHSTM(IMOL),LKINDM(IMOL)
          DO 4323 K=1,NRCM(IMOL)
            IF (ISRSM(IMOL).LT.3) THEN
              WRITE(IUNOT,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                          ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                          ISCDEM(IMOL,K),IESTMM(IMOL,K),
     .                          IBGKM(IMOL,K)
            ELSE IF (ISRSM(IMOL).EQ.3) THEN
              WRITE(IUNOT,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                          ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                          ISCD3M(IMOL,K),
     .                          ISCDEM(IMOL,K),IESTMM(IMOL,K),
     .                          IBGKM(IMOL,K)
            ELSE IF (ISRSM(IMOL).EQ.4) THEN
              WRITE(IUNOT,6666) IREACM(IMOL,K),IBULKM(IMOL,K),
     .                          ISCD1M(IMOL,K),ISCD2M(IMOL,K),
     .                          ISCD3M(IMOL,K),ISCD4M(IMOL,K),
     .                          ISCDEM(IMOL,K),IESTMM(IMOL,K),
     .                          IBGKM(IMOL,K)
            ENDIF
            WRITE(IUNOT,6664) EELECM(IMOL,K),EBULKM(IMOL,K),
     .                        ESCD1M(IMOL,K),ESCD2M(IMOL,K),
     .                        FREACM(IMOL,K)
 4323     CONTINUE
 4313   CONTINUE
*4C.   TEST IONS SPECIES CARDS, NIONI SPECIES
        WRITE(IUNOT,'(A72)') TXTION
        WRITE(IUNOT,6666) NIONI
        DO 4413 IION=1,NIONI
          IF (ISRFI(IION).GT.0 .AND. ISRFI(IION).LE.NATMI) THEN
            ISRFI(IION) = IATM1(ISRFI(IION))
          ELSE
            ISRFI(IION) = 0
          ENDIF
          IF (ABS(ISRTI(IION)).GT.0.AND.ABS(ISRTI(IION)).LE.NATMI) THEN
            ISRTI(IION) = SIGN(IATM1(ABS(ISRTI(IION))),ISRTI(IION))
          ELSE
            ISRTI(IION) = 0
          ENDIF
          IF (LKINDI(IION).GT.0 .AND. LKINDI(IION).LE.NATMI) THEN
            LKINDI(IION) = IATM1(LKINDI(IION))
          ELSE
            LKINDI(IION) = 0
          ENDIF
          WRITE(IUNOT,66666) IION,TEXTI(IION),NMASSI(IION),NCHARI(IION),
     .                         NPRTI(IION),NCHRGI(IION),
     .                         ISRFI(IION),ISRTI(IION),ISRSI(IION),
     .                         NRCI(IION),NFOLI(IION),NGENI(IION),
     .                         NHSTI(IION),LKINDI(IION)
          DO 4423 K=1,NRCI(IION)
            IF (ISRSI(IION).LT.3) THEN
              WRITE(IUNOT,6666) IREACI(IION,K),IBULKI(IION,K),
     .                          ISCD1I(IION,K),ISCD2I(IION,K),
     .                          ISCDEI(IION,K),IESTMI(IION,K),
     .                          IBGKI(IION,K)
            ELSE IF (ISRSI(IION).EQ.3) THEN
              WRITE(IUNOT,6666) IREACI(IION,K),IBULKI(IION,K),
     .                          ISCD1I(IION,K),ISCD2I(IION,K),
     .                          ISCD3I(IION,K),
     .                          ISCDEI(IION,K),IESTMI(IION,K),
     .                          IBGKI(IION,K)
            ELSE IF (ISRSI(IION).EQ.4) THEN
              WRITE(IUNOT,6666) IREACI(IION,K),IBULKI(IION,K),
     .                          ISCD1I(IION,K),ISCD2I(IION,K),
     .                          ISCD3I(IION,K),ISCD4I(IION,K),
     .                          ISCDEI(IION,K),IESTMI(IION,K),
     .                          IBGKI(IION,K)
            ENDIF
            WRITE(IUNOT,6664) EELECI(IION,K),EBULKI(IION,K),
     .                        ESCD1I(IION,K),ESCD2I(IION,K),
     .                        FREACI(IION,K)
 4423     CONTINUE
 4413   CONTINUE
*4D.   NEUTRAL PHOTONS SPECIES CARDS, NPHOTI SPECIES
        IF (TXTPHOT.NE.'' .OR. NPHOTI.GT.0) THEN
          WRITE(IUNOT,'(A72)') TXTPHOT
          WRITE(IUNOT,6666) NPHOTI
        END IF
        IF (NPHOTI.GT.0) THEN
          DO 4513 IPHOT=1,NPHOTI
            IF (ISRFPH(IPHOT).GT.0 .AND. ISRFPH(IPHOT).LE.NATMI) THEN
              ISRFPH(IPHOT) = IATM1(ISRFPH(IPHOT))
            ELSE
              ISRFPH(IPHOT) = 0
            ENDIF
            IF (ABS(ISRTPH(IPHOT)).GT.0.AND.
     &          ABS(ISRTPH(IPHOT)).LE.NATMI) THEN
              ISRTPH(IPHOT) =
     .         SIGN(IATM1(ABS(ISRTPH(IPHOT))),ISRTPH(IPHOT))
            ELSE
              ISRTPH(IPHOT) = 0
            ENDIF
            WRITE(IUNOT,66666) IPHOT,TEXTPH(IPHOT),0,0,0,0,
     .                         ISRFPH(IPHOT),ISRTPH(IPHOT),
     .                         ISRSPH(IPHOT),NRCPH(IPHOT),
     .                         NFOLPH(IPHOT),NGENPH(IPHOT),
     .                         NHSTPH(IPHOT)
            DO 4523 K=1,NRCPH(IPHOT)
              IF (ISRSPH(IPHOT).LT.3) THEN
                WRITE(IUNOT,6666) IREACPH(IPHOT,K),IBULKPH(IPHOT,K),
     .                            ISCD1PH(IPHOT,K),ISCD2PH(IPHOT,K),
     .                            ISCDEPH(IPHOT,K),IESTMPH(IPHOT,K),
     .                            IBGKPH(IPHOT,K)
              ELSE IF (ISRSPH(IPHOT).EQ.3) THEN
                WRITE(IUNOT,6666) IREACPH(IPHOT,K),IBULKPH(IPHOT,K),
     .                            ISCD1PH(IPHOT,K),ISCD2PH(IPHOT,K),
     .                            ISCD3PH(IPHOT,K),
     .                            ISCDEPH(IPHOT,K),IESTMPH(IPHOT,K),
     .                            IBGKPH(IPHOT,K)
              ELSE IF (ISRSPH(IPHOT).EQ.4) THEN
                WRITE(IUNOT,6666) IREACPH(IPHOT,K),IBULKPH(IPHOT,K),
     .                            ISCD1PH(IPHOT,K),ISCD2PH(IPHOT,K),
     .                            ISCD3PH(IPHOT,K),ISCD4PH(IPHOT,K),
     .                            ISCDEPH(IPHOT,K),IESTMPH(IPHOT,K),
     .                            IBGKPH(IPHOT,K)
              ENDIF
              WRITE(IUNOT,6664) EELECPH(IPHOT,K),EBULKPH(IPHOT,K),
     .                          ESCD1PH(IPHOT,K),ESCD2PH(IPHOT,K),
     .                          FREACPH(IPHOT,K),FLDLMPH(IPHOT,K)
 4523       CONTINUE
 4513     CONTINUE
        ENDIF
*4E. DATA FOR COLLISIONAL-RADIATIVE MODELS
        IF (NCRM.GT.0) THEN
          WRITE(IUNOT,'(A72)') TXTCRM
          WRITE(IUNOT,6666) NCRM
          DO ICRM=1,NCRM
            WRITE(IUNOT,6666) ICRM,NSCRM(ICRM),NIONREAC(ICRM),
     .                        NRCMBREAC(ICRM),TYPCRM(ICRM)
            DO K=1,NSCRM(ICRM)
              WRITE(IUNOT,6666) K,SCRM(K,ICRM),
     r                          NNCRM(K,ICRM),PHCRM(K,ICRM)
            END DO
          END DO
        ENDIF
*** 5. DATA FOR PLASMA BACKGROUND
        WRITE(IUNOT,'(A72)') TITLE5
*5A.   BULK ION SPECIES CARDS, NPLSI1 SPECIES
        WRITE(IUNOT,'(A72)') TXTPLS
        DO I = 1, IREAD-1
          WRITE(IUNOT,'(A72)') ZEILEN(IREAD)
        END DO
        WRITE(IUNOT,6666) NPLSI1
        DO 5114 IB=1,NPLSI1
          WRITE(IUNOT,66666) IB,TEXTP1(IB),NMASSP1(IB),NCHARP1(IB),
     .                                     NPRTP1(IB),NCHRGP1(IB),
     .                                     ISRFP1(IB),ISRTP1(IB),
     .                                     ISRSP1(IB),NRCP1(IB),
     .                                     DENSLIM1(IB),0,
     .                                     NHSTP1(IB),0,
     .                                     CDENMODEL1(IB),NRE1(IB)
          IF(ISP(IB).NE.0) THEN
            DO K=1, NRCP(ISP(IB))
              IF (ISRSP1(IB).LT.3) THEN
                WRITE(IUNOT,6666) IREACP(ISP(IB),K),IBULKP(ISP(IB),K),
     .                            ISCD1P(ISP(IB),K),ISCD2P(ISP(IB),K),
     .                            ISCDEP(ISP(IB),K)
              ELSE IF (ISRSP1(IB).EQ.3) THEN
                WRITE(IUNOT,6666) IREACP(ISP(IB),K),IBULKP(ISP(IB),K),
     .                            ISCD1P(ISP(IB),K),ISCD2P(ISP(IB),K),
     .                            ISCD3P(ISP(IB),K),ISCDEP(ISP(IB),K)
              ELSE IF (ISRSP1(IB).EQ.4) THEN
                WRITE(IUNOT,6666) IREACP(ISP(IB),K),IBULKP(ISP(IB),K),
     .                            ISCD1P(ISP(IB),K),ISCD2P(ISP(IB),K),
     .                            ISCD3P(ISP(IB),K),ISCD4P(ISP(IB),K),
     .                            ISCDEP(ISP(IB),K)
              ENDIF
              WRITE(IUNOT,6664) EELECP(ISP(IB),K),EBULKP(ISP(IB),K),
     .                          ESCD1P(ISP(IB),K),ESCD2P(ISP(IB),K),
     .                          FREACP(ISP(IB),K)
            END DO
          END IF
          SELECT CASE (CDENMODEL1(IB))
          CASE ('FORT.13   ')
            WRITE(IUNOT,6666) TDMISP1(1,IB)
          CASE ('BOLTZMANN ')
            WRITE(IUNOT,'(I6,6X,2E12.4)') TDMISP1(1,IB),
     .            TDMBOLTZ1(IB), TDMDELTAE1(IB)
          CASE ('CORONA    ')
            WRITE(IUNOT,'(I6,1X,A6,1X,A4,A9,A3,E12.4)')
     .            TDMISP1(1,IB),TDMFNAME(1,IB),
     .            TDMH21(1,IB),TDMREACTION(1,IB),
     .            TDMCR1(1,IB),TDM_CORONA(IB)
          CASE ('COLRAD    ')
            DO I = 1, NRE1(IB)
              WRITE(IUNOT,'(I6,1X,A6,1X,A4,A9,A3)')
     .              TDMISP1(I,IB),TDMFNAME(I,IB),
     .              TDMH21(I,IB),TDMREACTION(I,IB),
     .              TDMCR1(I,IB)
            END DO
          CASE DEFAULT
            CONTINUE
          END SELECT
 5114   CONTINUE
        deallocate(NMASSA,NMASSM,NMASSI,NMASSP,NMASSP1)
        deallocate(NCHARA,NCHARM,NCHARI,NCHARP,NCHARP1)
        deallocate(NPRTA,NPRTM,NPRTI,NPRTP,NPRTP1)
        deallocate(NCHRGI,NCHRGP,NCHRGP1)
        deallocate(ISRFA,ISRFM,ISRFI,ISRFP,ISRFP1)
        deallocate(ISRTA,ISRTM,ISRTI,ISRTP,ISRTP1)
        deallocate(ISRSA,ISRSM,ISRSI,ISRSP,ISRSP1)
        deallocate(NRCA,NRCM,NRCI,NRCP,NRCP1)
        deallocate(NFOLA,NFOLM,NFOLI)
        deallocate(NGENA,NGENM,NGENI)
        deallocate(NHSTA,NHSTM,NHSTI,NHSTP,NHSTP1)
        deallocate(TEXTA,TEXTM,TEXTI,TEXTP,TEXTP1)
        deallocate(IREACA,IBULKA,ISCD1A,ISCD2A,ISCD3A,ISCD4A,ISCDEA)
        deallocate(EELECA,EBULKA,ESCD1A,ESCD2A,FREACA,FLDLMA)
        deallocate(IREACM,IBULKM,ISCD1M,ISCD2M,ISCD3M,ISCD4M,ISCDEM)
        deallocate(EELECM,EBULKM,ESCD1M,ESCD2M,FREACM)
        deallocate(IREACI,IBULKI,ISCD1I,ISCD2I,ISCD3I,ISCD4I,ISCDEI)
        deallocate(EELECI,EBULKI,ESCD1I,ESCD2I,FREACI)
        deallocate(IREACP,IBULKP,ISCD1P,ISCD2P,ISCD3P,ISCD4P,ISCDEP)
        deallocate(EELECP,EBULKP,ESCD1P,ESCD2P,FREACP)
        deallocate(IESTMA,IBGKA,IESTMM,IBGKM,IESTMI,IBGKI)
        deallocate(DENSLIM,CDENMODEL,NRE,TDMISP)
        deallocate(DENSLIM1,CDENMODEL1,NRE1,TDMISP1)
        deallocate(TDMBOLTZ,TDMDELTAE,TDM_CORONA)
        deallocate(TDMBOLTZ1,TDMDELTAE1,TDM_CORONA1)
        deallocate(TDMFNAME,TDMH2,TDMREACTION,TDMCR)
        deallocate(TDMFNAME1,TDMH21,TDMREACTION1,TDMCR1)
        if (NPHOTI.GT.0) then
          deallocate(ISRFPH,ISRFP,ISRFP1)
          deallocate(ISRTPH,ISRTP,ISRTP1)
          deallocate(ISRSPH,ISRSP,ISRSP1)
          deallocate(NRCPH,NRCP,NRCP1)
          deallocate(NFOLPH)
          deallocate(NGENPH)
          deallocate(TEXTPH)
          deallocate(IREACPH,IBULKPH)
          deallocate(ISCD1PH,ISCD2PH,ISCD3PH,ISCD4PH,ISCDEPH)
          deallocate(EELECPH,EBULKPH,ESCD1PH,ESCD2PH,FREACPH,FLDLMPH)
          deallocate(IESTMPH,IBGKPH)
        endif
        if (NCRM.GT.0) then
          deallocate(NSCRM,NIONREAC,NRCMBREAC,TYPCRM)
          deallocate(SCRM,NNCRM,PHCRM)
        endif
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
            IF (IABS(INDPRO(2)).EQ.1 .OR. IABS(INDPRO(2)).EQ.3) THEN
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI0,PO,2)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI1,P1,2)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI2,P2,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI3,P3,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI4,P4,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI5,P5,3)
            ELSE IF (IABS(INDPRO(2)).EQ.2) THEN
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI0,PO,2)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI1,P1,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI2,P2,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI3,P3,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI4,P4,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI5,P5,3)
            ELSE
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI0,P0,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI1,P1,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI2,P2,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI3,P3,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI4,P4,3)
              call e_bundle_transmit_r
     .         (nplsi,nplsi1,e_bundle_contr,TI5,P5,3)
            ENDIF
            DO IB = 1, NPLSI1
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
          IF (IABS(INDPRO(3)).EQ.1 .OR. IABS(INDPRO(3)).EQ.3) THEN
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI0,PO,1)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI1,P1,1)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI5,P5,3)
          ELSE IF (IABS(INDPRO(3)).EQ.2) THEN
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI0,PO,1)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI1,P1,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI5,P5,3)
          ELSE
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI0,PO,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI1,P1,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,DI5,P5,3)
          ENDIF
          DO IB = 1, NPLSI1
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
          IF (IABS(INDPRO(4)).EQ.1 .OR. IABS(INDPRO(4)).EQ.3) THEN
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX0,PO,2)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX1,P1,2)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX5,P5,3)
          ELSE IF (IABS(INDPRO(4)).EQ.2) THEN
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX0,PO,2)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX1,P1,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX5,P5,3)
          ELSE
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX0,PO,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX1,P1,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VX5,P5,3)
          ENDIF
          DO IB = 1, NPLSI1
            WRITE(IUNOT,6664) P0(IB),P1(IB),P2(IB),
     .                        P3(IB),P4(IB),P5(IB)
          ENDDO
          deallocate(VX0,VX1,VX2,VX3,VX4,VX5)
          allocate(VY0(NPLSI),VY1(NPLSI),VY2(NPLSI),
     .             VY3(NPLSI),VY4(NPLSI),VY5(NPLSI))
          READ (IUNIN,6664) (VY0(I),VY1(I),VY2(I),VY3(I),VY4(I),VY5(I),
     .                       I=1,NPLSI)
          IF (IABS(INDPRO(4)).EQ.1 .OR. IABS(INDPRO(4)).EQ.3) THEN
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY0,PO,2)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY1,P1,2)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY5,P5,3)
          ELSE IF (IABS(INDPRO(4)).EQ.2) THEN
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY0,PO,2)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY1,P1,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY5,P5,3)
          ELSE
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY0,PO,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY1,P1,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VY5,P5,3)
          ENDIF
          DO IB = 1, NPLSI1
            WRITE(IUNOT,6664) P0(IB),P1(IB),P2(IB),
     .                        P3(IB),P4(IB),P5(IB)
          ENDDO
          deallocate(VY0,VY1,VY2,VY3,VY4,VY5)
          allocate(VZ0(NPLSI),VZ1(NPLSI),VZ2(NPLSI),
     .             VZ3(NPLSI),VZ4(NPLSI),VZ5(NPLSI))
          READ (IUNIN,6664) (VZ0(I),VZ1(I),VZ2(I),VZ3(I),VZ4(I),VZ5(I),
     .                       I=1,NPLSI)
          IF (IABS(INDPRO(4)).EQ.1 .OR. IABS(INDPRO(4)).EQ.3) THEN
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ0,PO,2)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ1,P1,2)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ5,P5,3)
          ELSE IF (IABS(INDPRO(4)).EQ.2) THEN
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ0,PO,2)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ1,P1,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ5,P5,3)
          ELSE
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ0,PO,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ1,P1,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ2,P2,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ3,P3,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ4,P4,3)
            call e_bundle_transmit_r
     .       (nplsi,nplsi1,e_bundle_contr,VZ5,P5,3)
          ENDIF
          DO IB = 1, NPLSI1
            WRITE(IUNOT,6664) P0(IB),P1(IB),P2(IB),
     .                        P3(IB),P4(IB),P5(IB)
          ENDDO
          deallocate(VZ0,VZ1,VZ2,VZ3,VZ4,VZ5)
        END IF
        deallocate(P0,P1,P2,P3,P4,P5)
        IF (INDPRO(5).EQ.5) THEN
          READ (IUNIN,'(A80)') ZEILE ! B0,B1,B2,B3,B4,B5
          WRITE(IUNOT,'(A80)') ZEILE
        ENDIF
        IF (INDPRO(12).EQ.5) THEN
          READ (IUNIN,'(A80)') ZEILE ! VL0,VL1,VL2,VL3,VL4,VL5
          WRITE(IUNOT,'(A80)') ZEILE
          IREAD=1
          IF (ZEILE(1:3) .NE. '***') IREAD = 0
        ENDIF
*** 6A. GENERAL DATA FOR REFLECTION MODEL
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        END IF
  610   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 610
        READ (ZEILE,6665) ldum(0) ! NLTRIM
        IREAD=0
        IF (ldum(0)) THEN
          READ (IUNIN,'(A80)') ZEILE
          IREAD=1
          IF (INDEX(ZEILE,'PATH')+INDEX(ZEILE,'path').EQ.0) THEN
615         IF (INDEX(ZEILE,'ON')+INDEX(ZEILE,'on').EQ.0) THEN
              GOTO 620
            ELSE
              WRITE(IUNOT,'(A80)') ZEILE
              READ (IUNIN,'(A80)') ZEILE
              IREAD=1
              GOTO 615
            ENDIF
          ELSE
            READ (ZEILE(7:lnblnk(ZEILE)),'(A74)') PATH
            I=INDEX(PATH,ACHAR(13)) !VK
            IF(I.GT.0.AND.I.LE.LEN_TRIM(PATH)) PATH(I:I)=' ' !VK TO REMOVE \r
            PATH=ADJUSTL(PATH)
            I2=INDEX(PATH,' ')
            FILE=REPEAT(' ',74)
            FILE(1:)=PATH(1:I2-1)
!            WRITE(IUNOT,'(A80)') ZEILE ! this writes the path
            READ (IUNIN,'(A80)') ZEILE
            IREAD=1
            NFR=0
  625       IF (INDEX(ZEILE,'ON')+INDEX(ZEILE,'on').EQ.0) THEN
              GOTO 620
            ELSE
              NFR=NFR+1
              READ (ZEILE,'(A74)') RFILNM
              I=INDEX(RFILNM,ACHAR(13)) !VK
              IF(I.GT.0.AND.I.LE.LEN_TRIM(RFILNM)) RFILNM(I:I)=' ' !VK TO REMOVE \r
              FILE(I2:)=RFILNM
              IF (NFR.EQ.1) THEN
                IF (open_file(29,FILE,
     .           'modules/Eirene/Database/Surfacedata/TRIM')) THEN  ! Path is valid, write it
                  WRITE(IUNOT,'(A6,A74)') 'cpath ',PATH
                  close(29)
                ELSE
                  WRITE(6,*) 
     .             'Invalid CPATH specification in this environment: '//
     .             'CPATH line removed'
                END IF
                WRITE(IUNOT,'(A74)') RFILNM
              ELSE
                WRITE(IUNOT,'(A74)') RFILNM
              END IF
              READ (IUNIN,'(A80)') ZEILE
              IREAD=1
              GOTO 625
            ENDIF
          ENDIF
        ENDIF
 620    allocate(DATD(NATMI),DATD1(NATMI1))
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,6664) (DATD(IATM),IATM=1,NATMI)
        ELSE IF (IREAD.EQ.1) THEN
          READ (ZEILE,6664) (DATD(IATM),IATM=1,MIN(NATMI,6))
          IREAD=0
          IF (NATMI.GT.6) THEN
            READ (IUNIN,6664) (DATD(IATM),IATM=7,NATMI)
          ENDIF
        ENDIF
        DATD1(1:NATMI1) = 0.0_R8
        DO IATM = 1, NATMI
          IF (IATM1(IATM).NE.0) THEN
            DATD1(IATM1(IATM)) = DATD(IATM)
          ENDIF
        ENDDO
        WRITE(IUNOT,6664) DATD1(1:NATMI1)
        deallocate(DATD,DATD1)
        READ (IUNIN,6664) rdum(1:NMOLI) ! (DMLD(IMOL),IMOL=1,NMOLI)
        WRITE(IUNOT,6664) rdum(1:NMOLI)
        READ (IUNIN,6664) rdum(1:NIONI) ! (DIOD(IION),IION=1,NIONI)
        WRITE(IUNOT,6664) rdum(1:NIONI)
        allocate(DPLD(NPLSI),DPLD1(NPLSI1))
        READ (IUNIN,6664) (DPLD(IPLS),IPLS=1,NPLSI)
        call e_bundle_transmit_r
     .   (nplsi,nplsi1,e_bundle_contr,DPLD,DPLD1,1)
        WRITE(IUNOT,6664) (DPLD1(IPLS),IPLS=1,NPLSI1)
        deallocate(DPLD,DPLD1)
        IF (NPHOTI.GT.0) THEN
          READ (IUNIN,6664) rdum(1:NPHOTI) ! (DPHD(IPHOT),IPHOT=1,NPHOTI)
          WRITE(IUNOT,6664) rdum(1:NPHOTI)
        END IF
        READ (IUNIN,'(A72)') ZEILE ! ERMIN,ERCUT,RPROB0,RINTEG,EINTEG,AINTEG
        WRITE(IUNOT,'(A72)') ZEILE
        DO
          IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A80)') ZEILE
            WRITE(IUNOT,'(A80)') ZEILE
          END IF
          IF (ZEILE(1:3) .EQ. '***') then
            IREAD = 1
            EXIT
          ENDIF
          IF (ZEILE(1:8) /= 'SURFMOD_') CYCLE
          IREAD=0
          READ (IUNIN,6666) idum(0:4) ! REFCUR%JLREF,REFCUR%JLSPT,
                                      ! REFCUR%JSRS,REFCUR%JSRC,
                                      ! REFCUR%JLCHSPNWL
          IF (idum(2).gt.0 .and. idum(2).le.NATMI)
     >      idum(2) = iatm1(idum(2))
          IF (idum(3).gt.0 .and. idum(3).le.NATMI) THEN
            idum(3) = iatm1(idum(3))
          ELSE IF (idum(3).gt.NATMI .and. idum(3).le.NATMI+NMOLI) THEN
            idum(3) = idum(3) - NATMI + NATMI1
          ENDIF
          IF (streql(eirene_format,'juelich') .or.
     .        streql(eirene_format,'new')) idum(4) = 1
          WRITE(IUNOT,6666) idum(0:4)
          READ (IUNIN,'(A80)') ZEILE ! REFCUR%ZNMLR,REFCUR%EWALLR,REFCUR%EWBINR,
                                     ! REFCUR%TRANSPR(1:2),REFCUR%FSHEATR
          WRITE(IUNOT,'(A80)') ZEILE
          READ (IUNIN,'(A80)') ZEILE ! REFCUR%RCYCFR,REFCUR%RCYCTR,REFCUR%RCPRMR,
                                     ! REFCUR%EXPPLR,REFCUR%EXPELR,REFCUR%EXPILR
          WRITE(IUNOT,'(A80)') ZEILE
          READ (IUNIN,'(A80)') ZEILE
          IREAD=1
          DO WHILE (ZEILE(1:1).NE.'*'.AND.ZEILE(1:8).NE.'SURFMOD_')
            varname = zeile(1:8)
            call EIRENE_uppercase (varname)
            select case (varname)
            case ('ISRS','ISRC','LCHSPNWL','TRANSP1','TRANSP2','RECYCF',
     .            'RECYCT','RECPRM','EXPPL','EXPEL','EXPIL','RECYCS',
     .            'RECYCC','SPTPRM','ESPUTS','ESPUTC')
              WRITE(IUNOT,'(A80)') ZEILE
              IREAD=0
            case ('ECHSPTG')
              WRITE(IUNOT,'(A8,A72)') 'ESPUTC  ',ZEILE(9:80)
              IREAD=0
            case ('FCHSPTG')
              IREAD=0
            case default
              IF (streql(eirene_format,'facelift')) THEN
                READ (ZEILE,6664) rdum(0:1), ! RECYCS(1,NLJ),RECYCC(1,NLJ),
     .                     rdum(4), rdum(5), ! ECHSPTG(1,NLJ),FCHSPTG(1,NLJ),
     .                     rdum(2)           ! SPTPRM(1,NLJ)
                rdum(3) = 0.0_R8
              ELSE IF (streql(eirene_format,'iter')) THEN
                READ (ZEILE,6664) rdum(0:4)  ! RECYCS(1,NLJ),RECYCC(1,NLJ),
                                             ! SPTPRM(1,NLJ),ESPUTS(1,NLJ),
                                             ! ESPUTC(1,NLJ)
              ELSE
                READ (ZEILE,6664) rdum(0:2)  ! RECYCS(1,NLJ),RECYCC(1,NLJ),
                                             ! SPTPRM(1,NLJ)
                rdum(3) = 0.0_R8
                rdum(4) = 0.0_R8
              ENDIF
              WRITE(IUNOT,6664) rdum(0:4)
              IREAD=0
            end select
            READ (IUNIN,'(A80)') ZEILE
            IREAD=1
          END DO
          IF (IREAD.EQ.1) WRITE(IUNOT,'(A80)') ZEILE
        END DO
*** 7. DATA FOR PRIMARY SOURCES, NSTRAI STRATA
cxpb  Here must read the whole block before rewrite
cxpb  to see if some strata need to be removed
  700   CONTINUE
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        ENDIF
  710   READ (IUNIN,'(A80)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') THEN
          WRITE(IUNOT,'(A80)') ZEILE
          GOTO 710
        END IF
        READ (ZEILE,6666) NSTRAI
        IREAD=0
        allocate (INDSRC(NSTRAI))
        allocate (TXTSOU(NSTRAI))
        allocate (NLAVRP(NSTRAI),NLAVRT(NSTRAI),
     .            NLSYMP(NSTRAI),NLSYMT(NSTRAI),
     .            NLRAY(NSTRAI))
        allocate (NPTS(NSTRAI),NINITL(NSTRAI),
     .            NEMODS(NSTRAI),NAMODS(NSTRAI),
     .            NMINPTS(NSTRAI),NPTSDEL(NSTRAI),
     .            NRAYEN(NSTRAI))
        allocate (FLUX(NSTRAI))
        allocate (NLATM(NSTRAI),NLMOL(NSTRAI),
     .            NLION(NSTRAI),NLPLS(NSTRAI),NLPHOT(NSTRAI))
        allocate (NSPEZ(NSTRAI))
        allocate (NLPNT(NSTRAI),NLLNE(NSTRAI),NLSRF(NSTRAI),
     .            NLVOL(NSTRAI), NLCNS(NSTRAI))
        allocate (NSRFSI(NSTRAI))
        allocate (INDIM(NSRFS,NSTRAI), INSOR(NSRFS,NSTRAI),
     .            INGRDA(NSRFS,NSTRAI,3),
     .            INGRDE(NSRFS,NSTRAI,3))
        allocate (SORWGT(NSRFS,NSTRAI), SORLIM(NSRFS,NSTRAI),
     .            SORIND(NSRFS,NSTRAI), SOREXP(NSRFS,NSTRAI),
     .            SORIFL(NSRFS,NSTRAI))
        allocate (NRSOR(NSRFS,NSTRAI), NPSOR(NSRFS,NSTRAI),
     .            NTSOR(NSRFS,NSTRAI), NBSOR(NSRFS,NSTRAI),
     .            NASOR(NSRFS,NSTRAI), NISOR(NSRFS,NSTRAI),
     .            ISTOR(NSRFS,NSTRAI))
        allocate (SORAD1(NSRFS,NSTRAI), SORAD2(NSRFS,NSTRAI),
     .            SORAD3(NSRFS,NSTRAI), SORAD4(NSRFS,NSTRAI),
     .            SORAD5(NSRFS,NSTRAI), SORAD6(NSRFS,NSTRAI))
        allocate (SORENI(NSTRAI), SORENE(NSTRAI),
     .            SORVDX(NSTRAI), SORVDY(NSTRAI), SORVDZ(NSTRAI))
        allocate (SORCOS(NSTRAI), SORMAX(NSTRAI),
     .            SORCTX(NSTRAI), SORCTY(NSTRAI), SORCTZ(NSTRAI))
        allocate (RAYFRAC(NSTRAI))
        allocate (SCALV(NSTRAI))
        allocate (IVLSF(NSTRAI),
     .            ISCLS(NSTRAI),ISCLT(NSTRAI),ISCL1(NSTRAI),
     .            ISCL2(NSTRAI),ISCL3(NSTRAI),ISCLB(NSTRAI),
     .            ISCLA(NSTRAI))
        savestratum(0) = .true.
        READ (IUNIN,6666) (INDSRC(ISTRA),ISTRA=1,NSTRAI)
        READ (IUNIN,6664) ALLOC, MPTS
        DO 712 ISTRA=1,NSTRAI
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
  713     READ (IUNIN,'(A80)') ZEILE
          IREAD=1
          IF (ZEILE(1:1) .EQ. '*') GOTO 713
          READ (ZEILE,6665)
     .     NLAVRP(ISTRA),NLAVRT(ISTRA),NLSYMP(ISTRA),NLSYMT(ISTRA),
     .     NLRAY(ISTRA)
          IREAD=0
          READ (IUNIN,6666) NPTS(ISTRA),NINITL(ISTRA),NEMODS(ISTRA),
     .                      NAMODS(ISTRA),NMINPTS(ISTRA),NPTSDEL(ISTRA),
     .                      NRAYEN(ISTRA)
          READ (IUNIN,66662) FLUX(ISTRA),SCALV(ISTRA),IVLSF(ISTRA),
     .                       ISCLS(ISTRA),ISCLT(ISTRA),ISCL1(ISTRA),
     .                       ISCL2(ISTRA),ISCL3(ISTRA),ISCLB(ISTRA),
     .                       ISCLA(ISTRA)
          ISCLS(ISTRA) = transform_tally(IVLSF(ISTRA).eq.1,ISCLT(ISTRA),
     .      ISCLS(ISTRA),iatm1,ipls1,nphoti,natmi,nmoli,nioni,nplsi,
     .                  natmi1,nplsi1)
          READ (IUNIN,6665) NLATM(ISTRA),NLMOL(ISTRA),
     .                      NLION(ISTRA),NLPLS(ISTRA),
     .                      NLPHOT(ISTRA)
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
            IF (streql(eirene_format,'old')) THEN
              IF (NLSRF(ISTRA)) THEN
                IF (INDIM(J,ISTRA).EQ.1.AND.
     &              INGRDA(J,ISTRA,1).NE.INGRDE(J,ISTRA,1)) THEN
                  INGRDA(J,ISTRA,2)=INGRDA(J,ISTRA,1)
                  INGRDE(J,ISTRA,2)=INGRDE(J,ISTRA,1)
                  INGRDA(J,ISTRA,1)=INSOR(J,ISTRA)
                  INGRDE(J,ISTRA,1)=INSOR(J,ISTRA)
                ENDIF
              ENDIF
            ENDIF
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
     .                      SORCTX(ISTRA),SORCTY(ISTRA),SORCTZ(ISTRA),
     .                      RAYFRAC(ISTRA)
  712   CONTINUE
cxpb  Start rewrite of block 7
        IS = 0
        do ISTRA=1,NSTRAI
          if (savestratum(istra)) then
            is = is + 1
            indsrc(is) = indsrc(istra)
          endif
        enddo
        NSTRA1 = is
        savestratum(NSTRAI+1:NSTRA) = .true.
        WRITE(IUNOT,6666) NSTRA1
        WRITE(IUNOT,6666) (INDSRC(ISTRA),ISTRA=1,NSTRA1)
        WRITE(IUNOT,6664) ALLOC, MPTS
        DO 7121 ISTRA=1,NSTRAI
          IF (INDSRC(ISTRA).EQ.6
     &     .or. .not.savestratum(ISTRA)) GOTO 7121
          WRITE(IUNOT,'(A72)') TXTSOU(ISTRA)
          WRITE(IUNOT,6665)
     .      NLAVRP(ISTRA),NLAVRT(ISTRA),NLSYMP(ISTRA),NLSYMT(ISTRA),
     .      NLRAY(ISTRA)
          WRITE(IUNOT,6666) NPTS(ISTRA),NINITL(ISTRA),NEMODS(ISTRA),
     .                      NAMODS(ISTRA),NMINPTS(ISTRA),NPTSDEL(ISTRA),
     .                      NRAYEN(ISTRA)
          WRITE(IUNOT,66662) FLUX(ISTRA),SCALV(ISTRA),IVLSF(ISTRA),
     .                       ISCLS(ISTRA),ISCLT(ISTRA),ISCL1(ISTRA),
     .                       ISCL2(ISTRA),ISCL3(ISTRA),ISCLB(ISTRA),
     .                       ISCLA(ISTRA)
          WRITE(IUNOT,6665) NLATM(ISTRA),NLMOL(ISTRA),
     .                      NLION(ISTRA),NLPLS(ISTRA),NLPHOT(ISTRA)
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
            IF (INDIM(J,ISTRA).EQ.1) THEN
              idum(3) = convert_eirene(non_commensurate,1,
     .           INSOR(J,ISTRA),(INGRDA(J,ISTRA,2)+INGRDE(J,ISTRA,2))/2,
     .           nx,ny,nx1,ny1,4)
              idum(4) = idum(3)
              idum(5) = idum(3)
              idum(6) = convert_eirene(non_commensurate,2,
     .           INSOR(J,ISTRA),INGRDA(J,ISTRA,2),nx,ny,nx1,ny1,1)
              idum(7) = convert_eirene(non_commensurate,2,
     .           INSOR(J,ISTRA),INGRDE(J,ISTRA,2),nx,ny,nx1,ny1,3)
            ELSE IF (INDIM(J,ISTRA).EQ.2) THEN
              idum(3) = convert_eirene(non_commensurate,2,
     .           (INGRDA(J,ISTRA,1)+INGRDE(J,ISTRA,1))/2,INSOR(J,ISTRA),
     .           nx,ny,nx1,ny1,4)
              idum(4) = convert_eirene(non_commensurate,1,
     .           INGRDA(J,ISTRA,1),INSOR(J,ISTRA),nx,ny,nx1,ny1,0)
              idum(5) = convert_eirene(non_commensurate,1,
     .           INGRDE(J,ISTRA,1),INSOR(J,ISTRA),nx,ny,nx1,ny1,2)
              idum(6) = idum(3)
              idum(7) = idum(3)
            ELSE IF (INDIM(J,ISTRA).EQ.3) THEN
              idum(4) = convert_eirene(non_commensurate,1,
     .           INGRDA(J,ISTRA,1),INGRDA(J,ISTRA,2),nx,ny,nx1,ny1,0)
              idum(5) = convert_eirene(non_commensurate,1,
     .           INGRDE(J,ISTRA,1),INGRDE(J,ISTRA,2),nx,ny,nx1,ny1,2)
              idum(6) = convert_eirene(non_commensurate,2,
     .           INGRDA(J,ISTRA,1),INGRDA(J,ISTRA,2),nx,ny,nx1,ny1,1)
              idum(7) = convert_eirene(non_commensurate,2,
     .           INGRDE(J,ISTRA,1),INGRDE(J,ISTRA,2),nx,ny,nx1,ny1,3)
            ELSE
              idum(4) = INGRDA(J,ISTRA,1)
              idum(5) = INGRDE(J,ISTRA,1)
              idum(6) = INGRDA(J,ISTRA,2)
              idum(7) = INGRDE(J,ISTRA,2)
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
     .                      SORCTX(ISTRA),SORCTY(ISTRA),SORCTZ(ISTRA),
     .                      RAYFRAC(ISTRA)
 7121   CONTINUE
        deallocate (INDSRC)
        deallocate (TXTSOU)
        deallocate (NLAVRP,NLAVRT,NLSYMP,NLSYMT,NLRAY)
        deallocate (NPTS,NINITL,NEMODS,NAMODS,NMINPTS,NPTSDEL,NRAYEN)
        deallocate (FLUX)
        deallocate (NLATM,NLMOL,NLION,NLPLS,NLPHOT)
        deallocate (NSPEZ)
        deallocate (NLPNT,NLLNE,NLSRF,NLVOL,NLCNS)
        deallocate (NSRFSI)
        deallocate (INDIM,INSOR,INGRDA,INGRDE)
        deallocate (SORWGT,SORLIM,SORIND,SOREXP,SORIFL)
        deallocate (NRSOR,NPSOR,NTSOR,NBSOR,NASOR,NISOR,ISTOR)
        deallocate (SORAD1,SORAD2,SORAD3,SORAD4,SORAD5,SORAD6)
        deallocate (SORENI,SORENE,SORVDX,SORVDY,SORVDZ)
        deallocate (SORCOS,SORMAX,SORCTX,SORCTY,SORCTZ)
        deallocate (RAYFRAC)
        deallocate (SCALV)
        deallocate (IVLSF,ISCLS,ISCLT,ISCL1,ISCL2,ISCL3,ISCLB,ISCLA)
*** 8. ADDITIONAL DATA FOR SPECIFIC ZONES
        READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
  810   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IREAD=1
        IF (ZEILE(1:1) .EQ. '*') GOTO 810
        READ (ZEILE,6666) NZADD
        IREAD=0
        DO 811 I=1,NZADD
          IF (IREAD.EQ.0) THEN
            READ (IUNIN,'(A80)') ZEILE
            WRITE(IUNOT,'(A80)') ZEILE
          ENDIF
          IREAD=0
          READ (IUNIN,6666) idum(0:1) ! INI,INE
          WRITE(IUNOT,6666) idum(0:1)
  812     READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
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
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        ENDIF
        IREAD=0
910     READ (IUNIN,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') THEN
          WRITE(IUNOT,'(A72)') ZEILE
          GOTO 910
        ELSE
          allocate(NLPRCA(NATM),NLPRCM(NMOL),NLPRCI(NION),
     .             NLPRCPH(NPHOT))
          allocate(NLPRCA1(NATM))
          READ (ZEILE,6665) (NLPRCA(J),J=1,NATM),(NLPRCM(J),J=1,NMOL),
     .                      (NLPRCI(J),J=1,NION),(NLPRCPH(J),J=1,NPHOT)
          NLPRCA1(1:NATM) = .false.
          DO IATM = 1, NATMI
            IF (IATM1(IATM).NE.0) THEN
              NLPRCA1(IATM1(IATM)) = NLPRCA(IATM)
            ENDIF
          ENDDO
          WRITE(IUNOT,6665) (NLPRCA1(J),J=1,NATM),(NLPRCM(J),J=1,NMOL),
     .                      (NLPRCI(J),J=1,NION),(NLPRCPH(J),J=1,NPHOT)
          deallocate(NLPRCA,NLPRCM,NLPRCI,NLPRCPH)
          deallocate(NLPRCA1)
        END IF
        IF (streql(eirene_format,'old')) THEN
          READ (IUNIN,'(A80)') ZEILE ! NLPRCF
          READ (IUNIN,'(A80)') ZEILE ! NLPRCD
          NPRCSF = 0
        ELSE
          READ (IUNIN,6666) NPRCSF
        ENDIF
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
C  DATA FOR BIAS SAMPLING
        READ (IUNIN,'(A72)') ZEILE ! WMINV,WMINS,WMINC,WMINL
        WRITE(IUNOT,'(A72)') ZEILE
        READ (IUNIN,'(A72)') ZEILE ! SPLPAR
        WRITE(IUNOT,'(A72)') ZEILE
C  DATA FOR STANDARD DEVIATION
        READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
cxpb Read first, transform and rewrite
        READ (IUNIN,6666) NSIGVI,NSIGSI,NSIGCI,
     .                    NSIGI_BGK,NSIGI_COP,NSIGI_SPC
        allocate (IGH(NSIGVI), IIH(NSIGVI))
        allocate (IGHW(NSIGSI), IIHW(NSIGSI))
        allocate (IGHC(2,NSIGCI), IIHC(2,NSIGCI))
        allocate (savesigv(NSIGVI))
        allocate (savesigs(NSIGSI))
        allocate (savesigc(NSIGCI))
        NSIGV1 = 0
        DO 913 J=1,NSIGVI
          READ (IUNIN,6666) IGH(J),IIH(J)
          IGH(J) = transform_tally(.true.,IIH(J),IGH(J),iatm1,ipls1,
     .                                 nphoti,natmi,nmoli,nioni,nplsi,
     .                                 natmi1,nplsi1)
          savesigv(J) = IGH(J).NE.0
          if (savesigv(J)) NSIGV1 = NSIGV1 + 1
  913   CONTINUE
        NSIGS1 = 0
        DO 914 J=1,NSIGSI
          READ (IUNIN,6666) IGHW(J),IIHW(J)
          IGHW(J) = transform_tally(.false.,IIHW(J),IGHW(J),
     .                              iatm1,ipls1,nphoti,natmi,nmoli,
     .                              nioni,nplsi,natmi1,nplsi1)
          savesigs(J) = IGHW(J).NE.0
          if (savesigs(J)) NSIGS1 = NSIGS1 + 1
  914   CONTINUE
        NSIGC1 = 0
        DO 915 J=1,NSIGCI
          READ (IUNIN,6666) IGHC(1,J),IIHC(1,J),IGHC(2,J),IIHC(2,J)
          IGHC(1,J) = transform_tally(.true.,IIHC(1,J),IGHC(1,J),
     .                                 iatm1,ipls1,nphoti,natmi,nmoli,
     .                                 nioni,nplsi,natmi1,nplsi1)
          IGHC(2,J) = transform_tally(.true.,IIHC(2,J),IGHC(2,J),
     .                                 iatm1,ipls1,nphoti,natmi,nmoli,
     .                                 nioni,nplsi,natmi1,nplsi1)
          savesigc(J) = IGHC(1,J).NE.0 .AND. IGHC(2,J).NE.0
          if (savesigc(J)) NSIGC1 = NSIGC1 + 1
  915   CONTINUE
        WRITE(IUNOT,6666) NSIGV1,NSIGS1,NSIGC1,
     .                    NSIGI_BGK,NSIGI_COP,NSIGI_SPC
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
** DATA FOR SOURCE CORRECTION
  920   READ (IUNIN,'(A80)') ZEILE
        IF (ZEILE(1:3) .EQ. '***') GOTO 1000
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 920
        K = NATM*NSTRA
        allocate(nlpbls(k),nlpbl1(k))
        IF (K.LE.60) THEN
          READ (ZEILE,6665) (nlpbls(i),i=1,k)
        ELSE
          KMOD = K/60
          DO J = 1, KMOD
            IF (j.eq.1) THEN
              READ (ZEILE,6665) (nlpbls(i),i=(j-1)*60+1,j*60)
            ELSE
              READ (IUNIN,6665) (nlpbls(i),i=(j-1)*60+1,j*60)
            ENDIF
          END DO
          READ (IUNIN,6665) (nlpbls(i),i=kmod*60+1,k)
        END IF
        ISTRA = 0
        nlpbl1(1:k) = .false.
        DO j = 1, nstrai
          IF (savestratum(j)) istra = istra+1
          DO i = 1, natmi
            IF (iatm1(i).NE.0 .AND. savestratum(j))
     .        nlpbl1(natm*(istra-1)+iatm1(i))=nlpbls(natm*(j-1)+i)
          END DO
        END DO
        WRITE(IUNOT,6665) nlpbl1(1:k)
        deallocate(nlpbls,nlpbl1)
*** 10. DATA FOR ADDITIONAL TALLIES, COLLISION ESTIMATORS AND ALGEBRAIC EXPRESSIONS
        READ (IUNIN,'(A80)') ZEILE
 1000   WRITE(IUNOT,'(A80)') ZEILE
 1010   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1010
        READ (ZEILE,6666) NADVI,NCLVI,NALVI,NADSI,NALSI,NADSPC
*** 10A. DATA FOR ADDITIONAL TALLIES
        READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        DO 1020 J=1,NADVI
 1021     READ (IUNIN,'(A80)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A80)') ZEILE
            GOTO 1021
          ELSE
            READ (ZEILE,6666) idum(0:3) ! IADVE(J),IADVS(J),IADVT(J),IADRC(J)
            idum(1) = transform_tally(.true.,idum(2),idum(1),
     .                                iatm1,ipls1,nphoti,natmi,nmoli,
     .                                nioni,nplsi,natmi1,nplsi1)
            WRITE(IUNOT,6666) idum(0:3)
          ENDIF
          READ (IUNIN,'(A72)') ZEILE ! TXTTLA(J)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
          WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 1020   CONTINUE
*** 10B. DATA FOR COLLISION ESTIMATORS
        READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        DO 1030 J=1,NCLVI
 1031     READ (IUNIN,'(A80)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A80)') ZEILE
            GOTO 1031
          ELSE
            READ (ZEILE,6666) idum(0:3) ! ICLVE(J),ICLVS(J),ICLVT(J),ICLRC(J)
            idum(1) = transform_tally(.true.,idum(2),idum(1),
     .                                iatm1,ipls1,nphoti,natmi,nmoli,
     .                                nioni,nplsi,natmi1,nplsi1)
            WRITE(IUNOT,6666) idum(0:3)
          ENDIF
          READ (IUNIN,'(A72)') ZEILE ! TXTTLC(J)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
          WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 1030   CONTINUE
*** 10C. DATA FOR ALGEBRAIC EXPRESSIONS
        READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        DO 1040 M=1,NALVI
 1041     READ (IUNIN,'(A80)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A80)') ZEILE
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
              idum(1) = transform_tally(.true.,idum(2),idum(1),
     .                                 iatm1,ipls1,nphoti,natmi,nmoli,
     .                                 nioni,nplsi,natmi1,nplsi1)
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
        READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        DO 1050 M=1,NADSI
 1051     READ (IUNIN,'(A80)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A80)') ZEILE
            GOTO 1051
          ELSE
            READ (ZEILE,6666) idum(0:3) ! IADSE(J),IADSS(J),IADST(J),IADSC(J)
            idum(1) = transform_tally(.false.,idum(2),idum(1),
     .                                iatm1,ipls1,nphoti,natmi,nmoli,
     .                                nioni,nplsi,natmi1,nplsi1)
            WRITE(IUNOT,6666) idum(0:3)
          ENDIF
          READ (IUNIN,'(A72)') ZEILE ! TXTTLW(J)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
          WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
 1050   CONTINUE
*** 10E. DATA FOR ALGEBRAIC SURFACE TALLIES
        READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        DO 1060 M=1,NALSI
 1061     READ (IUNIN,'(A80)') ZEILE
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A80)') ZEILE
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
              idum(1) = transform_tally(.false.,idum(2),idum(1),
     .                                  iatm1,ipls1,nphoti,natmi,nmoli,
     .                                  nioni,nplsi,natmi1,nplsi1)
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
*** 10F. DATA FOR SPECTRUM CALCULATION
        READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IREAD = 1
        IF (ZEILE(1:3).EQ.'** ' .AND. ZEILE(1:4).NE.'*** ') THEN
          DO M=1,NADSPC
            DO
              READ (IUNIN,'(A80)') ZEILE
              IF (ZEILE(1:1).EQ.'*') THEN
                WRITE(IUNOT,'(A80)') ZEILE
              ELSE
                EXIT
              ENDIF
            END DO
            READ (ZEILE,6666) idum(0:6) ! ISPSRF, IPTYP, IPSP, ISPTYP, NSPS,
                                            ! ISRFCLL, IDIREC
            IF (idum(1).EQ.1) idum(1) = iatm1(idum(1))
            IF (idum(1).EQ.4) idum(2) = ipls1(idum(2))
            WRITE(IUNOT,6666) idum(0:7)
            READ (IUNIN,'(A80)') ZEILE ! SPCMN, SPCMX, SPC_SHIFT,
                                       ! SPCPLT_X,SPCPLT_Y,SPCPLT_SAME
            WRITE(IUNOT,'(A80)') ZEILE
            IF (idum(6).NE.0) THEN
              READ (IUNIN,'(A80)') ZEILE ! SPCVX, SPCVY, SPCVZ
              WRITE(IUNOT,'(A80)') ZEILE
            END IF
          END DO
          IREAD = 0
        END IF
*** 11. DATA FOR NUMERICAL AND GRAPHICAL OUTPUT
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        END IF
 1110   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
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
        READ (IUNIN,6666) NVOLPR, idum(0) ! NSPCPR
        WRITE(IUNOT,6666) NVOLPR, idum(0)
        DO 1120 J=1,NVOLPR
          READ (IUNIN,6666) idum(0:4) ! NTLV,NFLGV,NSPZV1,NSPZV2,NTLVF
          idum(2) = transform_tally(.true.,idum(0),idum(2),iatm1,ipls1,
     .                                 nphoti,natmi,nmoli,nioni,nplsi,
     .                                 natmi1,nplsi1)
          idum(3) = transform_tally(.true.,idum(0),idum(3),iatm1,ipls1,
     .                                 nphoti,natmi,nmoli,nioni,nplsi,
     .                                 natmi1,nplsi1)
          WRITE(IUNOT,6666) idum(0:4)
 1120   CONTINUE
        READ (IUNIN,6666) NSURPR
        WRITE(IUNOT,6666) NSURPR
        DO 1130 J=1,NSURPR
          READ (IUNIN,6666) idum(0:4) ! NTLS,NFLGS,NSPZS1,NSPZS2,NTLSF
          idum(2) = transform_tally(.false.,idum(0),idum(2),iatm1,ipls1,
     .                                 nphoti,natmi,nmoli,nioni,nplsi,
     .                                 natmi1,nplsi1)
          idum(3) = transform_tally(.false.,idum(0),idum(3),iatm1,ipls1,
     .                                 nphoti,natmi,nmoli,nioni,nplsi,
     .                                 natmi1,nplsi1)
          WRITE(IUNOT,6666) idum(0:4)
 1130   CONTINUE
 1131   READ (IUNIN,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') THEN
          WRITE(IUNOT,'(A80)') ZEILE
          GOTO 1131
        ELSE
          READ (ZEILE,6665) ldum(0:14) ! PL1ST,PL2ND,PL3RD,PLADD,PLHST,
                                       ! PLCUT(1),PLCUT(2),PLCUT(3),PLBOX,PLSTOR,
                                       ! PLNUMV,PLNUMS,PLARR,LRPSCUT,PLIDL
          LRPSCUT = ldum(13)
          IF (streql(eirene_format,'old')) THEN
            ldum(11) = ldum(10)        ! PLNUMB -> PLNUMS
            ldum(10) = .false.         ! no PLNUMV in "old" Eirene
          ENDIF
          WRITE (IUNOT,6665) ldum(0:14)
        ENDIF
        READ (IUNIN,6666) NPLINR,NPLOTR,NPLDLR,NPLINP,NPLOTP,NPLDLP,
     .                    NPLINT,NPLOTT,NPLDLT,idum(9) ! nsrfcls
        IF (NPLOTR.LE.0) NPLOTR=NR1ST
        IF (NPLOTP.LE.0) NPLOTP=NP2ND
        IF (NPLOTT.LE.0) NPLOTT=NT3RD
        NPLOTR=MIN0(NPLOTR,NR1ST)
        if (nltri.or.nltet) NPLOTR=MIN0(NPLOTR,NR1ST-1)
        NPLINR=MAX0(NPLINR,1)
        NPLDLR=MAX0(NPLDLR,1)
        NPLOTP=MIN0(NPLOTP,NP2ND)
        NPLINP=MAX0(NPLINP,1)
        NPLDLP=MAX0(NPLDLP,1)
        NPLOTT=MIN0(NPLOTT,NT3RD)
        NPLINT=MAX0(NPLINT,1)
        NPLDLT=MAX0(NPLDLT,1)
        idum(0) = convert_eirene
     .   (non_commensurate,1,NPLINR,NPLINP,nx,ny,nx1,ny1,0)
        idum(1) = convert_eirene
     .   (non_commensurate,1,NPLOTR,NPLOTP,nx,ny,nx1,ny1,2)
        IF (NPLDLR.ne.1) THEN
          idum(2) = max(1, (NPLDLR*(ny1-1))/ny)
        ELSE ! we assume here that all new cells are to be plotted
          idum(2) = 1
        ENDIF
        idum(3) = convert_eirene
     .   (non_commensurate,2,NPLINR,NPLINP,nx,ny,nx1,ny1,1)
        idum(4) = convert_eirene
     .   (non_commensurate,2,NPLOTR,NPLOTP,nx,ny,nx1,ny1,3)
        IF (NPLDLP.ne.1) THEN
          idum(5) = max(1,
     .     (NPLDLP*(nx1-max(nncut*3,1)))/(nx-max(nncut*3,1)))
        ELSE ! we assume here that all new cells are to be plotted
          idum(5) = 1
        ENDIF
        WRITE(IUNOT,6666) idum(0:5),NPLINT,NPLOTT,NPLDLT,idum(9)
        do j=1,abs(idum(9)) !{
          read (iunin,6666) idum(0)  ! msrfcls(j)
          if(idum(0).gt.0) then !{
            read (iunin,6666) (idum(i),i=1,idum(0))
          end if !}
        end do !}
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
        DO I=2,9 ! Convert to new symbol numbers
          IF (streql(eirene_format,'old')) THEN
            IF (idum(I).eq.4) THEN
              idum(I) = 6
            ELSE IF (idum(I).ge.6 .and. idum(I).le.8) THEN
              idum(I) = idum(I)+2
            ELSE IF (idum(I).ge.10 .and. idum(I).le.14) THEN
              idum(I) = idum(I)+2
            ELSE IF (idum(I).eq.15) THEN
              idum(I) = 18
            END IF
          ELSE IF (streql(eirene_format,'new')) THEN
            IF (idum(I).ge.4 .and. idum(I).le.14) THEN
              idum(I) = idum(I)+2
            ELSE IF (idum(I).eq.15) THEN
              idum(I) = 18
            END IF
          ELSE IF (streql(eirene_format,'facelift')) THEN
            IF (idum(I).ge.3 .and. idum(I).le.16) THEN
              idum(I) = idum(I)+2
            ENDIF
          ENDIF
        END DO
        WRITE(IUNOT,6666) idum(0:10)
C  DATA FOR PLOTS OF VOLUME TALLIES
 1151   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:1) .EQ. '*') GOTO 1151
        READ (ZEILE,6666) NVOLPL
        IF (NVOLPL.LE.0) GOTO 1175
        READ (IUNIN,6665) (ldum(J),J=0,NSTRA)
        WRITE(IUNOT,6665) (ldum(J),J=0,NSTRA)
        IF (streql(eirene_format,'iter') .or.
     .      streql(eirene_format,'juelich')) then
          IF (LRPSCUT) THEN
            READ (IUNIN,'(A72)') ZEILE ! CUTPLANE(1:4)
            WRITE(IUNOT,'(A72)') ZEILE
          END IF
        ENDIF
        DO 1150 J=1,NVOLPL
          IF (streql(eirene_format,'iter') .or.
     .        streql(eirene_format,'juelich') .or.
     .        streql(eirene_format,'facelift') .or.
     .        streql(eirene_format,'new')) then
 1152       READ (IUNIN,'(A80)') ZEILE
            WRITE(IUNOT,'(A80)') ZEILE
            IF (ZEILE(1:1) .EQ. '*') GOTO 1152
            READ (ZEILE,6666) NSP
          ELSE IF (streql(eirene_format,'old')) then
            READ (IUNIN,6666) NTL,NSP
            WRITE(IUNOT,6666) NSP
          ENDIF
          READ (IUNIN,6665) ldum(0:3) ! PLTL2D(J),PLTL3D(J),PLTLLG(J),PLTLER(J)
          WRITE(IUNOT,6665) ldum(0:3)
          READ (IUNIN,'(A72)') ZEILE ! TALZMI(J),TALZMA(J),TALXMI(J),TALXMA(J),
                                     ! TALYMI(J),TALYMA(J)
          WRITE(IUNOT,'(A72)') ZEILE
          IF (ldum(0)) THEN
            READ (IUNIN,6665) ldum(2:3) ! LHIST2(J),LSMOT2(J)
            WRITE(IUNOT,6665) ldum(2:3)
            DO 1160 I=1,NSP
              IF (streql(eirene_format,'iter') .or.
     .            streql(eirene_format,'juelich') .or.
     .            streql(eirene_format,'facelift') .or.
     .            streql(eirene_format,'new')) then
                READ (IUNIN,6666) ISPL,NTL,NPLIN2,NPLOT2,NPLDL2
              ELSE IF (streql(eirene_format,'old')) then
                READ (IUNIN,6666) ISPL,NPLIN2,NPLOT2,NPLDL2
              ENDIF
              idum(1) = transform_tally(.true.,NTL,ISPL,iatm1,ipls1,
     .                                 nphoti,natmi,nmoli,nioni,nplsi,
     .                                 natmi1,nplsi1)
              idum(2) = NTL
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
            IF (streql(eirene_format,'iter') .or.
     .          streql(eirene_format,'juelich') .or.
     .          streql(eirene_format,'facelift') .or.
     .          streql(eirene_format,'new')) then
              READ (IUNIN,6665) ldum(2:10) ! LHIST3(J),LCNTR3(J),LSMOT3(J),
                                           ! LRAPS3(J),LVECT3(J),LRPVC3(J),
                                           ! LRPS3D,LRPSCN
              WRITE(IUNOT,6665) ldum(2:10)
              READ (IUNIN,6665) ldum(2:4) ! LPRAD3(J),LPPOL3(J),LPTOR3(J)
              WRITE(IUNOT,6665) ldum(2:4)
            ELSE IF (streql(eirene_format,'old')) then
              READ (IUNIN,6665) ldum(2:8) ! LHIST3(J),LCTOUR(J),LSMOT3(J),
                                          ! LRAPS(J),LPRAD(J),LPPOL(J),LPTOR(J)
              WRITE(IUNOT,6665) ldum(2:5),
     .                          .false.,.false.,.false.,.false.
              WRITE(IUNOT,6665) ldum(6:8)
            ENDIF
            DO 1161 I=1,NSP
              IF (streql(eirene_format,'iter') .or.
     .            streql(eirene_format,'juelich') .or.
     .            streql(eirene_format,'facelift') .or.
     .            streql(eirene_format,'new')) then
                READ (IUNIN,6666)
     &           ISPL,NTL,IPROJ,NPLI1,NPLO1,NPLI2,NPLO2,idum(8)  ! IPLN
              ELSE IF (streql(eirene_format,'old')) then
                READ (IUNIN,6666) ISPL,IPROJ,NPLI1,NPLO1,NPLI2,NPLO2
                idum(8) = 0
              ENDIF
              NPLI1 = MAX(1,NPLI1)
              NPLO1 = MIN(NR1ST,NPLO1)
              NPLI2 = MAX(1,NPLI2)
              NPLO2 = MIN(NP2ND,NPLO2)
              idum(1) = transform_tally(.true.,NTL,ISPL,iatm1,ipls1,
     .                                 nphoti,natmi,nmoli,nioni,nplsi,
     .                                 natmi1,nplsi1)
              idum(2) = NTL
              IF (ldum(2)) THEN  ! LPRAD
                idum(3) = convert_eirene(non_commensurate,
     .           1,IPROJ,(NPLI2+NPLO2)/2,nx,ny,nx1,ny1,4)
              ELSE IF (ldum(3)) THEN  ! LPPOL
                idum(3) = convert_eirene(non_commensurate,
     .           2,(NPLI1+NPLO1)/2,IPROJ,nx,ny,nx1,ny1,4)
              ELSE IF (ldum(4)) THEN  ! LPTOR
                idum(3) = IPROJ
              ENDIF
              idum(4) = convert_eirene(non_commensurate,
     .         1,NPLI1,NPLI2,nx,ny,nx1,ny1,0)
              idum(5) = convert_eirene(non_commensurate,
     .         1,NPLO1,NPLO2,nx,ny,nx1,ny1,2)
              idum(6) = convert_eirene(non_commensurate,
     .         2,NPLI1,NPLI2,nx,ny,nx1,ny1,1)
              idum(7) = convert_eirene(non_commensurate,
     .         2,NPLO1,NPLO2,nx,ny,nx1,ny1,3)
              WRITE(IUNOT,6666) idum(1:8)
 1161       CONTINUE
            READ (IUNIN,'(A72)') ZEILE ! TALW1(J),TALW2(J),FCABS1(J),FCABS2(J),RPSDL
            WRITE(IUNOT,'(A72)') ZEILE
          ENDIF
 1150   CONTINUE
 1175   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IREAD=1
        IF (ZEILE(1:3) .EQ. '***') GOTO 1205
        GOTO 1175
*** 12. DATA FOR DIAGNOSTIC MODULE
 1205   CONTINUE
 1210   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IREAD=0
        IF (ZEILE(1:1) .EQ. '*') GOTO 1210
        READ (ZEILE,6666) NCHORI,idum(0),lntspec,lntsput,idum(1),
     .                    nlstadds,idum(2),ngwties
        IF (NCHORI.LE.0) GOTO 1230
        DO 1220 I=1,NCHORI
          READ (IUNIN,'(A72)') ZEILE ! TXTSIG(I)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,6666) idum(0:2) ! NCHTAL(I),NSPSCL(I),NSPNEW(I)
          WRITE(IUNOT,6666) idum(0:2)
          READ (IUNIN,6666) idum(0:6) ! NSPSTR(I),NSPSPZ(I),
                                      ! NSPINI(I),NSPEND(I),
                                      ! NSPBLC(I),NSPADD(I),MASLOC
          if (.not.savestratum(idum(0))) idum(0) = 0
          idum(1) = iatm1(idum(1))
          WRITE(IUNOT,6666) idum(0:6)
          READ (IUNIN,'(A72)') ZEILE ! EMIN1(I),EMAX1(I)
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(A72)') ZEILE ! IPIVOT(I),XPIVOT,YPIVOT,ZPIVOT
          WRITE(IUNOT,'(A72)') ZEILE
          READ (IUNIN,'(A72)') ZEILE ! ICHORD(I),XCHORD,YCHORD,ZCHORD
          WRITE(IUNOT,'(A72)') ZEILE
 1220   CONTINUE
        READ (IUNIN,6665) ldum(0:1) ! PLCHOR,PLSPEC
        WRITE(IUNOT,6665) ldum(0:1)
 1225   CONTINUE
        IF (lntspec.GT.0) THEN
          READ (IUNIN,'(2e12.0,4i6)') rdum(0:1), idum(0:3) ! emnspc,emxspc,
                                                           ! npspcen,iatspec,npspctt,npspcfi
          WRITE(IUNOT,'(2e12.0,4i6)') rdum(0:1), idum(0:3)
        END IF
        DO I=1,LNTSPUT !{
          READ (IUNIN,'(A80)') ZEILE  ! listmtmp(i),txtmtmp(i),
                                      ! listmtra(i),txtmtra(i)
          WRITE(IUNOT,'(A80)') ZEILE
        END DO !}
        IF (NLSTADDS.GT.0) THEN !{
          READ (IUNIN,6666) (idum(i),i=1,nlstadds) ! (listadds(i),i=1,nlstadds)
          WRITE(IUNOT,6666) (idum(i),i=1,nlstadds)
          READ (IUNIN,6666) (idum(i),i=1,nlstadds) ! (listaddg(i),i=1,nlstadds)
          WRITE(IUNOT,6666) (idum(i),i=1,nlstadds)
        END IF !}
        DO J=1,NGWTIES !{
          READ (IUNIN,6666) idum(1), idum(0) ! kgwties(j),mgwties(j)
          WRITE(IUNOT,6666) idum(1), idum(0)
          IF (idum(1).GT.0) THEN !{
            READ (IUNIN,6666) (idum(i),i=1,idum(0)) ! (lgwties(i,j),i=1,mgwties(j))
            WRITE(IUNOT,6666) (idum(i),i=1,idum(0))
          END IF !}
        END DO !}
 1227   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:3).EQ.'CH4') GOTO 1227
        IREAD=1
 1230   CONTINUE
CVK DATA FOR POST-PROCESSING WINDOWS
        READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        IF (ZEILE(1:3).NE.'***') THEN !{
          NPWIN=0
          READ(ZEILE,6666) NPWIN
          IF (NPWIN.GT.0) THEN  !{
            DO I=1,NPWIN
              IF (IREAD.GT.0) THEN
                IREAD=1
             ELSE
               READ (IUNIN,6666) idum(1),idum(0) ! J,NPWSURF(I)
               WRITE(IUNOT,6666) idum(1),idum(0)
             END IF
             READ (IUNIN,6666) (idum(k),k=1,idum(0)) ! (SURF_BUF(ISURF+K),K=1,NPWSURF(I))
             WRITE(IUNOT,6666) (idum(k),k=1,idum(0))
             READ (IUNIN,6666) idum(1),idum(0) ! J,NPWTALLIE(I)
             WRITE(IUNOT,6666) idum(1),idum(0)
             READ (IUNIN,6666) (idum(k),k=1,idum(0)) ! (TALLIES_BUF(ITALLLIE+K),K=1,NPWTALLIE(I))
             DO K = 1,idum(0)
               KMOD = MOD(idum(K),100) ! PWSPECIE
               KINT = INT(idum(K)/100) ! PWTALLIE
               KMOD = transform_tally(.true.,KINT,KMOD,iatm1,ipls1,
     .                                 nphoti,natmi,nmoli,nioni,nplsi,
     .                                 natmi1,nplsi1)
               idum(K) = kint*100 + kmod
             END DO
             WRITE(IUNOT,6666) (idum(k),k=1,idum(0))
             READ (IUNIN,'(A80)') ZEILE
             WRITE(IUNOT,'(A80)') ZEILE
             IREAD=1
             IF (ZEILE(1:3).EQ.'///') THEN
               READ (IUNIN,6666) (idum(k),k=1,idum(0)) ! (TALLIES_D_BUF(ITALLLIE+K),K=1,NPWTALLIE(I))
               DO K = 1,idum(0)
                 KMOD = MOD(idum(K),100) ! PWSPECIE_D
                 KINT = INT(idum(K)/100) ! PWTALLIE_D
                 KMOD = transform_tally(.false.,KINT,KMOD,iatm1,ipls1,
     .                                  nphoti,natmi,nmoli,nioni,nplsi,
     .                                  natmi1,nplsi1)
                 idum(K) = kint*100 + kmod
               END DO
               WRITE(IUNOT,6666) (idum(k),k=1,idum(0))
             END IF
            END DO
          END IF !}
        END IF !}
        GOTO 1232
C  SKIP READING REST OF THIS BLOCK
 1231   READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
 1232   IREAD=0
        IF (ZEILE(1:3).NE.'***') GOTO 1231
        IREAD=1
*** 13. DATA FOR ITERATIVE AND TIME DEP. OPTION
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        END IF
        IREAD=0
        READ (IUNIN,6666) idum(0) ! NPRNLI
        WRITE(IUNOT,6666) idum(0)
        IF (idum(0).LE.0) GOTO 1350
        READ (IUNIN,'(A80)') ZEILE ! NPTST,NTMSTP
        WRITE(IUNOT,'(A80)') ZEILE
        IREAD=1
        IF (ZEILE(1:1).EQ.'*') GOTO 1350
        IREAD=0
        READ (IUNIN,'(A80)') ZEILE ! DTIMV,TIME0
        WRITE(IUNOT,'(A80)') ZEILE
*** 13A. DATA FOR SNAPSHOT TALLIES
        READ (IUNIN,'(A80)') ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        READ (IUNIN,6666) NSNVI
        WRITE(IUNOT,6666) NSNVI
        DO 1320 J=1,NSNVI
 1321     READ (IUNIN,'(A80)') ZEILE
          IREAD=1
          IF (ZEILE(1:1) .EQ. '*') THEN
            WRITE(IUNOT,'(A80)') ZEILE
            GOTO 1321
          ELSE
            READ (ZEILE,6666) idum(0:3) ! ISNVE(J),ISNVS(J),ISNVT(J),ISNRC(J)
            IREAD=0
            idum(1) = transform_tally(.true.,idum(2),idum(1),
     .                                iatm1,ipls1,nphoti,natmi,nmoli,
     .                                nioni,nplsi,natmi1,nplsi1)
            WRITE(IUNOT,6666) idum(0:3)
            READ (IUNIN,'(A72)') ZEILE ! TXTTAL(J,NTALT)
            WRITE(IUNOT,'(A72)') ZEILE
            READ (IUNIN,'(2A24)') TXTSPC,TXTUNT
            WRITE(IUNOT,'(2A24)') TXTSPC,TXTUNT
          END IF
 1320   CONTINUE
C  SKIP READING REST OF THIS BLOCK
 1350   CONTINUE
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        END IF
        IF (ZEILE(1:3).NE.'***') GOTO 1350
        IREAD=1
*** 14. DATA FOR INTERFACING ROUTINE "INFCOP"
        IF (IREAD.EQ.0) THEN
          READ (IUNIN,'(A80)') ZEILE
          WRITE(IUNOT,'(A80)') ZEILE
        END IF
        IF (NMODE.EQ.0) THEN
          READ (IUNIN,6666) NAINI,idum(1:2) ! NCOPII,NCOPIE
          WRITE(IUNOT,6666) NAINI,idum(1:2)
          DO 3020 J=1,NAINI
 3021       READ (IUNIN,'(A80)') ZEILE
            WRITE(IUNOT,'(A80)') ZEILE
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
        ELSE IF (NMODE.NE.0) THEN
          READ (IUNIN,'(A72)') ZEILE ! LSYMET,LBALAN,LCHKQUD
          WRITE(IUNOT,'(A72)') ZEILE
          IF (streql(eirene_format,'old')) READ (IUNIN,'(A72)') ZEILE
          READ (IUNIN,'(9I6)') idum(0:8) ! NFLA,NCUTB,NCUTL,IMF,
                                         ! ntrfrm,nfull,ibrad,ibpol,ibtor
          idum(0) = ns1-natmi1
          if (idum(6).eq.0) idum(6)=1
          if (idum(7).eq.0) idum(7)=1
          if (idum(8).eq.0) idum(8)=1
          WRITE(IUNOT,'(9I6)') idum(0:8)
          allocate(IFLB1(NPLSI1),FCTE1(NPLSI1),BMASS1(NPLSI1),
     .             LKINDP1(NPLSI1))
          DO I=1,NPLSI
            READ (IUNIN,'(2I6,2E12.4,I6)')
     .       idum(0),IFLB(I),FCTE(I),BMASS(I),LKINDP(I)
          ENDDO
          IPLS = 0
          DO IATM=1,NATMI1
            IS = EB2SPCR1(IATM)+1
            IF (IATM.EQ.NATMI1) THEN
              JS = NS1-1
            ELSE
              JS = EB2SPCR1(IATM+1)-1
            ENDIF
            DO I = IS, JS
              IPLS = IPLS + 1
              IFLB1(IPLS) = IPLS
              LKINDP1(IPLS) = IATM
            ENDDO
          ENDDO
          call e_bundle_transmit_r
     .     (nplsi,nplsi1,e_bundle_contr,FCTE,FCTE1,3)
          call e_bundle_transmit_r
     .     (nplsi,nplsi1,e_bundle_contr,BMASS,BMASS1,3)
          DO I=1,NPLSI
            IF (I.GT.(NS-NATMI)) THEN
              IFLB1(I-(NS-NATMI)+(NS1-NATMI1))=IFLB(I)
              FCTE1(I-(NS-NATMI)+(NS1-NATMI1))=FCTE(I)
              BMASS1(I-(NS-NATMI)+(NS1-NATMI1))=BMASS(I)
              LKINDP1(I-(NS-NATMI)+(NS1-NATMI1))=LKINDP(I)
            ELSEIF (ipls1(I).GT.0) THEN
              IF (LKINDP(I).GT.0 .AND. LKINDP(I).LE.NATMI) THEN
                LKINDP1(IPLS1(I))=IATM1(LKINDP(I))
              END IF
            END IF
          END DO
          DO I=1,NPLSI1
            WRITE(IUNOT,'(1P,2I6,2E12.4,I6)')
     .       I,IFLB1(I),FCTE1(I),BMASS1(I),LKINDP1(I)
          END DO
          deallocate(IFLB1,FCTE1,BMASS1,LKINDP1)
          READ (IUNIN,'(2I6)') idum(0:1) ! NDXA,NDYA
          WRITE(IUNOT,'(2I6)') nx1+max(1,nncut*3)-1, ny1
C  NUMBER OF TARGETS ON B2 SURFACES: NTARGI
          READ (IUNIN,'(I6)') NTARGI
          WRITE(IUNOT,'(I6)') NTARGI
          allocate(NTGPRT(NTARGI))
          IF (NTARGI.GT.0) THEN
            IF (streql(eirene_format,'old')) THEN
              NTGPRT(1:NTARGI) = 1
            ELSE
              READ (IUNIN,6666) (NTGPRT(I),I=1,NTARGI)
            ENDIF
            WRITE(IUNOT,6666) (NTGPRT(I),I=1,NTARGI)
          END IF
          DO 3030 I=1,NTARGI
            DO 3033 J=1,NTGPRT(I)
 3031         READ (IUNIN,'(A80)') ZEILE
              IREAD=1
              IF (ZEILE(1:1).EQ.'*') THEN
                WRITE(IUNOT,'(A80)') ZEILE
                GOTO 3031
              ENDIF
              READ (ZEILE,6666) idum(0),NDT(I,J),NINCT(I,J),NIXY(I,J),
     .                          NTIN(I,J),NTEN(I,J),NIFLG(I,J),
     .                          NPTC(I,J),NPTCM(I,J),
     .                          NSPZI(I,J),NSPZE(I,J),idum(10)
              NSPZI(I,J)=MAX0(1,NSPZI(I,J))
              NSPZE(I,J)=MIN0(NS-NATMI,NSPZE(I,J))
              IF (NSPZE(I,J).LT.NSPZI(I,J)) THEN
                NSPZI(I,J)=1
                NSPZE(I,J)=NS-NATMI
              ENDIF
              idum(8) = is1beginrange(NSPZI(I,J),NSPZE(I,J),ns,ns1,is1)
              idum(9) = is1endrange(NSPZI(I,J),NSPZE(I,J),ns,ns1,is1)
              IF (NIXY(I,J).EQ.1) THEN  ! target runs radially
                IF (NINCT(I,J).EQ.1) THEN
                  iside = 3
                ELSE IF (NINCT(I,J).EQ.-1) THEN
                  iside = 1
                ENDIF
                idum(1) = convert_eirene(non_commensurate,2,
     .           (NTIN(I,J)+NTEN(I,J))/2,NDT(I,J)+1,
     .           nx,ny,nx1,ny1,iside)-1
                idum(4) = convert_eirene(non_commensurate,1,
     .           NTIN(I,J),NDT(I,J)+1,nx,ny,nx1,ny1,0)
                idum(5) = convert_eirene(non_commensurate,1,
     .           NTEN(I,J),NDT(I,J)+1,nx,ny,nx1,ny1,2)
              ELSE IF (NIXY(I,J).EQ.2) THEN  ! target runs poloidally
                IF (NINCT(I,J).EQ.1) THEN
                  iside = 2
                ELSE IF (NINCT(I,J).EQ.-1) THEN
                  iside = 0
                ENDIF
                idum(1) = convert_eirene(non_commensurate,1,NDT(I,J)+1,
     .           (NTIN(I,J)+NTEN(I,J))/2,nx,ny,nx1,ny1,iside)-1
                idum(4) = convert_eirene(non_commensurate,2,NDT(I,J)+1,
     .           NTIN(I,J),nx,ny,nx1,ny1,1)
                idum(5) = convert_eirene(non_commensurate,2,NDT(I,J)+1,
     .           NTEN(I,J),nx,ny,nx1,ny1,3)
              ENDIF
              WRITE(IUNOT,6666) idum(0:1),NINCT(I,J),NIXY(I,J),
     .         idum(4:5),NIFLG(I,J),NPTC(I,J),NPTCM(I,J),idum(8:10)
 3033       CONTINUE
 3030     CONTINUE
          READ (IUNIN,'(A80)') ZEILE ! CHGP,CHGEE,CHGEI,CHGMOM
          WRITE(IUNOT,'(A80)') ZEILE
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
            idum(1) = transform_tally(.false.,idum(2),idum(1),
     .                                iatm1,ipls1,nphoti,natmi,nmoli,
     .                                nioni,nplsi,natmi1,nplsi1)
            WRITE(IUNOT,'(6I6)') idum(0:2)
 3050     CONTINUE
        ENDIF
*** 15. ADDITIONAL GEOMETRY
        IF (NMODE.GT.0) THEN
          READ (IUNIN,'(A80)',end=8686) ZEILE  ! geometry_comment
          WRITE(IUNOT,'(A80)') ZEILE
          READ (IUNIN,'(A80)') ZEILE
          DO WHILE (ZEILE(1:1).EQ.'*') !{
            WRITE(IUNOT,'(A80)') ZEILE
            READ (IUNIN,'(A80)') ZEILE
          END DO !}
          IF (streql(eirene_format,'juelich')) THEN    ! using geousr.F from user_general
            READ (ZEILE,6666) NADMOD,NASMOD,NORMOD
            WRITE(IUNOT,6666) NADMOD,NASMOD,NORMOD
            IGEOUSR = 0
          ELSE IF (streql(eirene_format,'iter')) THEN  ! using geousr.F from user_iter
            IF (INDEX(ZEILE,'GENERAL') /= 0) THEN      ! we read only nadmod, nasmod, normod
              WRITE(IUNOT,'(A80)') ZEILE
              READ (IUNIN,6666) NADMOD,NASMOD,NORMOD
              WRITE(IUNOT,6666) NADMOD,NASMOD,NORMOD
              IGEOUSR = 0
            ELSEIF (INDEX(ZEILE,'BIASED_GARCHING') /= 0) THEN ! we read only onetwo, limpos (reverse order !)
              WRITE(IUNOT,'(A80)') ZEILE
              NADMOD = 0
              NASMOD = max(npplg/3,1)*4
              NORMOD = 0
              IGEOUSR = -1
              WRITE(IUNOT,6666) NADMOD,NASMOD,NORMOD
              WRITE(IUNOT,'(A80)') ZEILE
            ELSE             ! we read both, but limpos, onetwo ! This will be the default writing format
              READ (ZEILE,6666) NADMOD,NASMOD,NORMOD
              WRITE(IUNOT,6666) NADMOD,NASMOD,NORMOD
              IGEOUSR = 1
            END IF
          ELSE IF (streql(eirene_format,'facelift')) THEN  ! we read both, but limpos, onetwo, xpt, ypt
            READ (ZEILE,6666) NADMOD,NASMOD,NORMOD
            WRITE(IUNOT,6666) NADMOD,NASMOD,NORMOD
            IGEOUSR = 1
          ELSE IF (streql(eirene_format,'old') .or.
     .             streql(eirene_format,'new')) THEN ! we read only onetwo, limpos (reverse order !)
            NADMOD = 0
            NASMOD = max(npplg/3,1)*4
            NORMOD = 0
            IGEOUSR = -1
            WRITE(IUNOT,6666) NADMOD,NASMOD,NORMOD
          END IF

          DO I=1,NADMOD
            READ (IUNIN,'(A80)') ZEILE ! NRS,IPUNKT,XCOOR,YCOOR,ZCOOR
            WRITE(IUNOT,'(A80)') ZEILE
          ENDDO
          DO I=1,NASMOD
            IF (streql(eirene_format,'old') .or.
     .          streql(eirene_format,'new')) THEN
              IF (I.EQ.1) THEN
                READ (ZEILE,6666) idum(0:1)  ! onetwo(i),limpos(i)
                WRITE(IUNOT,6666) idum(1), idum(0)  ! limpos(i), onetwo(i)
              ELSE
                READ (IUNIN,6666) idum(0:1)  ! onetwo(i),limpos(i)
                WRITE(IUNOT,6666) idum(1), idum(0)  ! limpos(i), onetwo(i)
              END IF
              IF (idum(0).lt.0) THEN
                READ (IUNIN,6666) idum(2:3)   ! xpolpos(i), ypolpos(i)
                WRITE(IUNOT,6666) idum(2:3)
              END IF
            ELSE IF (streql(eirene_format,'facelift')) THEN
              READ (IUNIN,66665) idum(0:1), rdum(0:1) ! limpos(i), onetwo(i), xpt, ypt
              WRITE(IUNOT,6666) idum(0:1)
            ELSE IF (streql(eirene_format,'juelich')) THEN
              READ (IUNIN,6666) idum(3:6) ! NAS,IPUNKT,NSSIR,NSSIP
              NSSIR = idum(5)
              NSSIP = idum(6)
              idum(5) = convert_eirene(non_commensurate,1,
     .         NSSIR,NSSIP,nx,ny,nx1,ny1,4)
              idum(6) = convert_eirene(non_commensurate,2,
     .         NSSIR,NSSIP,nx,ny,nx1,ny1,4)
              WRITE(IUNOT,6666) idum(3:6)
            ELSE IF (streql(eirene_format,'iter')) THEN
              IF (IGEOUSR == 0) THEN  ! 'GENERAL' format
                READ (IUNIN,6666) idum(3:6) ! NAS,IPUNKT,NSSIR,NSSIP
                NSSIR = idum(5)
                NSSIP = idum(6)
                idum(5) = convert_eirene(non_commensurate,1,
     .           NSSIR,NSSIP,nx,ny,nx1,ny1,4)
                idum(6) = convert_eirene(non_commensurate,2,
     .           NSSIR,NSSIP,nx,ny,nx1,ny1,4)
                WRITE(IUNOT,6666) idum(3:6)
              ELSE IF (IGEOUSR == -1) THEN  ! 'BIASED_GARCHING' format
                READ (IUNIN,6666) idum(0:1)  ! onetwo(i),limpos(i)
                WRITE(IUNOT,6666) idum(0:1)  ! onetwo(i),limpos(i)
                IF (idum(0).lt.0) THEN
                  READ (IUNIN,6666) idum(2:3)   ! xpolpos(i), ypolpos(i)
                  WRITE(IUNOT,6666) idum(2:3)
                END IF
              ELSE IF (IGEOUSR == 1) THEN   ! default format
                READ (IUNIN,6666) idum(0:1)  ! limpos(i), onetwo(i)
                WRITE(IUNOT,6666) idum(0:1)  ! limpos(i), onetwo(i)
                IF (idum(1).lt.0) THEN
                  READ (IUNIN,6666) idum(2:3)   ! xpolpos(i), ypolpos(i)
                  WRITE(IUNOT,6666) idum(2:3)
                END IF
              END IF
            END IF
          ENDDO
          DO I=1,NORMOD
            READ (IUNIN,6666) idum(7:9) ! IDIR,IR,IP
            IR = idum(8)
            IP = idum(9)
            idum(8) = convert_eirene(non_commensurate,1,
     .       IR,IP,nx,ny,nx1,ny1,4)
            idum(9) = convert_eirene(non_commensurate,2,
     .       IR,IP,nx,ny,nx1,ny1,4)
            WRITE(IUNOT,6666) idum(7:9)
          ENDDO
        END IF
*** IMFORMATION FOR MPI
 8685   CONTINUE
        READ (IUNIN,'(A80)',end=8686) ZEILE
        WRITE(IUNOT,'(A80)') ZEILE
        goto 8685
 8686   continue
        close(IUNOT)
        goto 8688
 8687   continue
        write(*,*) 'Empty input.dat file not converted'
        savestratum(0:nstrai) = .true.
 8688   continue
        close(IUNIN)
cxpb  Convert fort.44 file if needed
        inquire(file='fort.44',exist=todo)
        if (todo.and..not.streql(solps_version,'iter')) then
          allocate(epls_map(ns1-natmi1,ns-natmi))
          epls_map(1:ns1-natmi1,1:ns-natmi)=
     &     e_bundle_contr(1:ns1-natmi1,1:ns-natmi)
          call b2yt_ngread(solps_version, eirene_format,
     &     nx1, ny1, natmi1, iatm1, ns1-natmi1, ipls1, ns-natmi,
     &     savestratum, non_commensurate, epls_map,
     &     jsep1, region1, leftcut1, rightcut1, bottomcut1, topcut1,
     &     bottomiy1,
     &     hx, hy, qc, hx1, hy1, qc1)
          write(nout(5),'(a)') 'fort.44.2'
          deallocate(epls_map)
        endif
      else
        iatm1(1:natm) = 0
        savestratum(0:nstrai) = .true.
! must fill new arrays b2espcr1 and eb2spcr1 not provided by modified input.dat file
        iss = 0
        do is = 1, nspecies
          if (is1(eb2spcr(is)).ne.-1) then
            iss = iss + 1
            eb2spcr1(iss) = is1(eb2spcr(is))
            b2espcr1(
     .   is1beginrange(eb2spcr(is),eb2spcr(is)+nfluids(is),ns,ns1,is1):
     .     is1endrange(eb2spcr(is),eb2spcr(is)+nfluids(is),ns,ns1,is1))=
     .       is
          end if
        end do
      endif
#else
      iatm1(1:natm) = 0
      savestratum(0:nstrai) = .true.
! must fill new arrays b2espcr1 and eb2spcr1 not provided by modified input.dat file
      iss = 0
      do is = 1, nspecies
        if (is1(eb2spcr(is)).ne.-1) then
          iss = iss + 1
          eb2spcr1(iss) = is1(eb2spcr(is))
          b2espcr1(
     .   is1beginrange(eb2spcr(is),eb2spcr(is)+nfluids(is),ns,ns1,is1):
     .     is1endrange(eb2spcr(is),eb2spcr(is)+nfluids(is),ns,ns1,is1))=
     .     is
        end if
      end do
#endif
! fill up the mapping between sequence numbers
      if (same_bundling) then
        do is = 1, nspecies
          ispecies1(is) = is
        end do
      else
        do is = 1, nspecies
          if (is1(eb2spcr(is)).ne.-1) then
            ispecies1(is) = b2espcr1(is1(eb2spcr(is)))
          endif
        end do
      end if

*  ..rewrite neutrals namelist
*  ..recall it was already read earlier
      inquire(file='../b2.neutrals.parameters',exist=todo)
      new_filename = '../b2.neutrals.parameters.2'
      inquire(file=new_filename,exist=found)
      if (found) write(6,*) trim(new_filename), ' already present!'
      do while((neutrals_namelist.eq.1.or.use_eirene.ne.0.or.todo)
     & .and..not.found)
        ib1 = 0
        nstrai1 = 0
        nstrat1 = 0
        allocate(istra1(nstrat))
        do ib = 1, nstrat
         if (savestratum(ib) .or. crcstra(ib).eq.'T') then
          ib1 = ib1 + 1
          istra1(ib) = ib1
          nstrat1 = nstrat1 + 1
          if (crcstra(ib).ne.'T') nstrai1 = nstrai1 + 1
         else
          istra1(ib) = 0
         end if
        end do
        allocate(recyc1(0:ns1-1,nstrat1))
        allocate(erecyc1(0:ns1-1,nstrat1))
        allocate(rcion1(0:ns1-1,nstrat1))
        allocate(phys_sput1(0:ns1-1,nstrat1))
        allocate(chem_sput1(0:ns1-1,nstrat1))
        allocate(lstrascl1(nstrat1,natm))
        do ib = 1, nstrat
         if (savestratum(ib) .or. crcstra(ib).eq.'T') then
          ib1 = istra1(ib)
          crcstra1(ib1) = crcstra(ib)
          if (crcstra(ib).eq.'S') then
            rcpos1(ib1) = conformal_convert(non_commensurate,
     .       2,rcstart(ib),rcpos(ib),nx,ny,nx1,ny1,2)
            rcstart1(ib1) = conformal_convert(non_commensurate,
     .       1,rcstart(ib),rcpos(ib),nx,ny,nx1,ny1,1)
            rcend1(ib1) = conformal_convert(non_commensurate,
     .       1,rcend(ib),rcpos(ib),nx,ny,nx1,ny1,3)
            if (rcpos1(ib).eq.-2) then
              rc_list_size1(ib) = 0
              do is = 1, rc_list_size(ib)
                ix1 = conformal_convert(non_commensurate,
     .           1,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,1)
                ix2 = conformal_convert(non_commensurate,
     .           1,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,3)
                iy1 = conformal_convert(non_commensurate,
     .           2,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,2)
                do i = ix1, ix2
                  j = 0
                  found = .false.
                  do while (j.le.rc_list_size1(ib) .and. .not.found)
                    j = j + 1
                    found = i.eq.rc_list_x1(j,ib) .and.
     &                    iy1.eq.rc_list_y1(j,ib)
                  end do
                  if (.not.found) then
                    rc_list_size1(ib) = rc_list_size1(ib) + 1
                    rc_list_x1(rc_list_size1(ib),ib) = i
                    rc_list_y1(rc_list_size1(ib),ib) = iy1
                  end if
                end do
              end do
            end if
          else if (crcstra(ib).eq.'N') then
            rcpos1(ib1) = conformal_convert(non_commensurate,
     .       2,rcstart(ib),rcpos(ib),nx,ny,nx1,ny1,0)
            rcstart1(ib1) = conformal_convert(non_commensurate,
     .       1,rcstart(ib),rcpos(ib),nx,ny,nx1,ny1,1)
            rcend1(ib1) = conformal_convert(non_commensurate,
     .       1,rcend(ib),rcpos(ib),nx,ny,nx1,ny1,3)
            if (rcpos1(ib).eq.-2) then
              rc_list_size1(ib) = 0
              do is = 1, rc_list_size(ib)
                ix1 = conformal_convert(non_commensurate,
     .           1,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,1)
                ix2 = conformal_convert(non_commensurate,
     .           1,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,3)
                iy1 = conformal_convert(non_commensurate,
     .           2,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,0)
                do i = ix1, ix2
                  j = 0
                  found = .false.
                  do while (j.le.rc_list_size1(ib) .and. .not.found)
                    j = j + 1
                    found = i.eq.rc_list_x1(j,ib) .and.
     &                    iy1.eq.rc_list_y1(j,ib)
                  end do
                  if (.not.found) then
                    rc_list_size1(ib) = rc_list_size1(ib) + 1
                    rc_list_x1(rc_list_size1(ib),ib) = i
                    rc_list_y1(rc_list_size1(ib),ib) = iy1
                  end if
                end do
              end do
            end if
          else if (crcstra(ib).eq.'W') then
            rcpos1(ib1) = conformal_convert(non_commensurate,
     .       1,rcpos(ib),rcstart(ib),nx,ny,nx1,ny1,3)
            rcstart1(ib1) = conformal_convert(non_commensurate,
     .       2,rcpos(ib),rcstart(ib),nx,ny,nx1,ny1,0)
            rcend1(ib1) = conformal_convert(non_commensurate,
     .       2,rcpos(ib),rcend(ib),nx,ny,nx1,ny1,2)
            if (rcpos1(ib).eq.-2) then
              rc_list_size1(ib) = 0
              do is = 1, rc_list_size(ib)
                ix1 = conformal_convert(non_commensurate,
     .           1,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,3)
                iy1 = conformal_convert(non_commensurate,
     .           2,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,0)
                iy2 = conformal_convert(non_commensurate,
     .           2,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,2)
                do j = iy1, iy2
                  i = 0
                  found = .false.
                  do while (i.le.rc_list_size1(ib) .and. .not.found)
                    i = i + 1
                    found = j.eq.rc_list_y1(j,ib) .and.
     &                    ix1.eq.rc_list_x1(j,ib)
                  end do
                  if (.not.found) then
                    rc_list_size1(ib) = rc_list_size1(ib) + 1
                    rc_list_x1(rc_list_size1(ib),ib) = ix1
                    rc_list_y1(rc_list_size1(ib),ib) = j
                  end if
                end do
              end do
            end if
          else if (crcstra(ib).eq.'E') then
            rcpos1(ib1) = conformal_convert(non_commensurate,
     .       1,rcpos(ib),rcstart(ib),nx,ny,nx1,ny1,1)
            rcstart1(ib1) = conformal_convert(non_commensurate,
     .       2,rcpos(ib),rcstart(ib),nx,ny,nx1,ny1,0)
            rcend1(ib1) = conformal_convert(non_commensurate,
     .       2,rcpos(ib),rcend(ib),nx,ny,nx1,ny1,2)
            if (rcpos1(ib).eq.-2) then
              rc_list_size1(ib) = 0
              do is = 1, rc_list_size(ib)
                ix1 = conformal_convert(non_commensurate,
     .           1,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,1)
                iy1 = conformal_convert(non_commensurate,
     .           2,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,0)
                iy2 = conformal_convert(non_commensurate,
     .           2,rc_list_x(is,ib),rc_list_y(is,ib),nx,ny,nx1,ny1,2)
                do j = iy1, iy2
                  i = 0
                  found = .false.
                  do while (i.le.rc_list_size1(ib) .and. .not.found)
                    i = i + 1
                    found = j.eq.rc_list_y1(j,ib) .and.
     &                    ix1.eq.rc_list_x1(j,ib)
                  end do
                  if (.not.found) then
                    rc_list_size1(ib) = rc_list_size1(ib) + 1
                    rc_list_x1(rc_list_size1(ib),ib) = ix1
                    rc_list_y1(rc_list_size1(ib),ib) = j
                  end if
                end do
              end do
            end if
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
          if (l_chem_sput) then
          cftmp(0:ns-1) = chem_sput(0:ns-1,ib)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          chem_sput1(0:ns1-1,ib1)=cftmp(0:ns1-1)
          chemsp1(ib1)=chemsp(ib)
          else
            chem_sput1(0:ns1-1,ib1) = 0.0_R8
            chemsp1(ib1) = .false.
          endif
          recyceir1(ib1)=recyceir(ib)
          volrecstart1(ib1)=volrecstart(ib)
          userfluxparm1(ib1,1:2)=userfluxparm(ib,1:2)
         endif
        enddo
cxpb We modify the values in stored in track_index before re-bundling
cxpb as they no longer necessarily point to the same species
        do is = 0, ns-1
          if (track_index(is).eq.1) then
            track_index(is) = itrack1
          elseif (track_index(is).eq.2) then
            track_index(is) = itrack2
          elseif (track_index(is).eq.3) then
            track_index(is) = itrack3
          else
            track_index(is) = 0
          endif
        end do
        call convert_bundle_i
     .   (ns1,ns,nsm,is_neutral1,bundle_frac,track_index)
        do is = 1, nspecies1
          if (eb2spcr1(is).ne.-1) then
            isp1(is) = b2espcr1(eb2spcr1(is))
          else
            isp1(is) = 0
          end if
        end do
        mlcmp1 = 0.0_R8
        do is = 1, nspecies1
          if (isp1(is).ne.0) then
            mlcmp1(isp1(is),1:nmol) = mlcmp(is,1:nmol)
          endif
        end do
        do is = 1, nmol
          if (lmolscl(is).ne.0) then
            if (iatm1(lmolscl(is)).gt.0 .and.
     .          iatm1(lmolscl(is)).le.nspecies)
     >                lmolscl(is) = ispecies1(lmolscl(is))
          end if
        end do
        do is = 1, nion
          if (lionscl(is).ne.0) then
            if (iatm1(lionscl(is)).gt.0 .and.
     .          iatm1(lionscl(is)).le.nspecies)
     >                lionscl(is) = ispecies1(lionscl(is))
          end if
        end do
        targsp1 = 0
        do is = 1, ntrack
          ib1 = 0
          do ib = 1, nstrat
            if (savestratum(ib).or.crcstra(ib).eq.'T') then
              ib1 = ib1+1
              if (targsp(ib,is).ge.0) targsp1(ib1,is)=is1(targsp(ib,is))
            endif
          enddo
        enddo
        write(*,*) 'Writing ',trim(new_filename)
        write(nout(5),'(a)') trim(new_filename)
        open(34,file=new_filename)
        write(34,'(a)') '&NEUTRALS'
        write(34,'(a,i2,a)') 'NSTRAI=', nstrai1,','
        if (eirene_step_cpu.gt.0.0_R8) write(34,'(a,1pg10.3,a)')
     >   'EIRENE_STEP_CPU=',eirene_step_cpu,','
        write(34,'(a,1pg10.3,a)')
     >   'EIRENE_STEP_DT= ',eirene_step_dt,','
        if (eirene_mod.gt.1)
     >   write(34,'(a,i2,a)') 'EIRENE_MOD=', eirene_mod,','
        line = 'CRCSTRA=      '
        len=14
        do ib = 1, nstrat1
          line = line(1:len)//'  '''//crcstra1(ib)//''','
          len=len+3+1+2
        end do
        write(34,'(a)') line(1:len)
        write(hlp_frm,'(a,i2,a)') '(a,',nstrat1,'(1x,i4,'',''))'
        if(minval(rcpos1(1:nstrat1)).gt.-2) then
          write(34, hlp_frm) 'RCPOS=        ',
     .     (rcpos1(ib)  ,ib=1,nstrat1)
          write(34, hlp_frm) 'RCSTART=      ',
     .     (rcstart1(ib),ib=1,nstrat1)
          write(34, hlp_frm) 'RCEND=        ',
     .     (rcend1(ib)  ,ib=1,nstrat1)
        else
          do ib = 1, nstrat1
            if(rcpos1(ib).gt.-2) then
              write(35,'(1x,a,i2,a,i4)') 'RCPOS(',ib,')=  ',rcpos1(ib)
              write(35,'(1x,a,i2,a,i4)') 'RCSTART(',ib,')=',rcstart1(ib)
              write(35,'(1x,a,i2,a,i4)') 'RCEND(',ib,')=  ',rcend1(ib)
            else
              write(35,'(1x,a,i2,a,i4)') 'RC_LIST_SIZE(',ib,')=',
     .                                    rc_list_size1(ib)
              write(hlp_frm,'(a,i4,a)')
     .         '(1x,a,i2,a,',rc_list_size1(ib),'(1x,i4,'',''))'
              write(35, hlp_frm) 'RC_LIST_X(1,',ib,')=',
     .         rc_list_x1(1:bc_list_size1(ib),ib)
              write(35, hlp_frm) 'RC_LIST_Y(1,',ib,')=',
     .         rc_list_y1(1:bc_list_size1(ib),ib)
            end if
          end do
        endif
        write(hlp_frm,'(a,i2,a)') '(a,',nstrat1,'(1x,i4,'',''))'
        if (maxval(species_start1(1:nstrat1)).gt.0)
     >   write(34, hlp_frm)
     .             'SPECIES_START=',(species_start1(ib),ib=1,nstrat1)
        if (minval(species_end1(1:nstrat1)).lt.ns1-1)
     >   write(34, hlp_frm)
     .             'SPECIES_END=  ',(species_end1(ib),ib=1,nstrat1)
        write(hlp_frm,'(a,i3,a)') '(a,i2,a,1p,',ns1,'(1x,g10.3,'',''))'
        do ib = 1, nstrat1
          write(34, hlp_frm)
     .     'RECYC(0,',ib,')= ',(recyc1(is,ib),is=0,ns1-1)
          write(34, hlp_frm)
     .     'ERECYC(0,',ib,')=',(erecyc1(is,ib),is=0,ns1-1)
          if (maxval(rcion1(0:ns1-1,ib)).gt.0.0_R8)
     >     write(34, hlp_frm)
     .     'RCION(0,',ib,')= ',(rcion1(is,ib),is=0,ns1-1)
        end do
        write(hlp_frm,'(a,i3,a)') '(a,i2,a,1p,',ns1,'(1x,g10.3,'',''))'
        todo = .false.
        do ib = 1, nstrat1
          todo = todo .or. chemsp1(ib)
        enddo
        write(hlp_frm,'(a,i2,a)') '(a,',nstrat1,'l1,a)'
        if (todo) write(34,hlp_frm) 'CHEMSP=',chemsp1(1:nstrat1),','
        write(hlp_frm,'(a,i3,a)') '(a,',ns1,'(1x,i2,'',''))'
        write(34, hlp_frm) 'TRACK_INDEX=',(track_index(is),is=0,ns1-1)
        track_chem_sput1 = .false.
        if (itrack1.gt.0) track_chem_sput1(itrack1) = track_chem_sput(1)
        if (itrack2.gt.0) track_chem_sput1(itrack2) = track_chem_sput(2)
        if (itrack3.gt.0) track_chem_sput1(itrack3) = track_chem_sput(3)
        write(hlp_frm,'(a,i2,a)') '(a,',itrack,'l1,a)'
        if (itrack.gt.0) write(34, hlp_frm)
     .   'TRACK_CHEM_SPUT=',track_chem_sput1(1:itrack),','
        write(hlp_frm,'(a,i2,a)') '(a,i2,a,',nstrat1,'(1x,i3,'',''))'
        do is = 1, ntrack
          if (maxval(targsp1(1:nstrat1,is)).gt.0) write(34, hlp_frm)
     .     'TARGSP(1,',is,')=',(targsp1(ib,is),ib=1,nstrat1)
        end do
        write(hlp_frm,'(a,i2,a)') '(a,1p,',nstrat1,'(1x,g10.3,'',''))'
        if (minval(recyceir1(1:nstrat1)).lt.1.0_R8)
     >   write(34, hlp_frm) 'RECYCEIR=',(recyceir1(1:nstrat1))
        write(hlp_frm,'(a,i3,a)') '(a,i2,a,1p,',ns1,'(1x,g10.3,'',''))'
        do ib = 1, nstrat1
          if (maxval(phys_sput1(0:ns1-1,ib)).gt.0.0_R8)
     >     write(34, hlp_frm)
     .     'PHYS_SPUT(0,',ib,')=',(phys_sput1(is,ib),is=0,ns1-1)
          if (maxval(chem_sput1(0:ns1-1,ib)).gt.0.0_R8)
     >     write(34, hlp_frm)
     .     'CHEM_SPUT(0,',ib,')=',(chem_sput1(is,ib),is=0,ns1-1)
        end do
        if (chemical_erosion_redep_fac.ne.1.0_R8)
     >   write(34, '(a,1pg10.3,a)')
     .    'CHEMICAL_EROSION_REDEP_FAC=',
     .     chemical_erosion_redep_fac,','
        if (chemical_erosion_be_fac) then
          write(34, '(a)')
     >     'CHEMICAL_EROSION_BE_FAC=T,'
          if (chemical_erosion_be_fac_a.ne.0.0_R8)
     >     write(34, '(a,1pg10.3,a)')
     .     'CHEMICAL_EROSION_BE_FAC_A=',
     .      chemical_erosion_be_fac_a,','
          if (chemical_erosion_be_fac_b.ne.0.0_R8)
     >     write(34, '(a,1pg10.3,a)')
     .     'CHEMICAL_EROSION_BE_FAC_B=',
     .      chemical_erosion_be_fac_b,','
          if (chemical_erosion_be_fac_c.ne.0.0_R8)
     >     write(34, '(a,1pg10.3,a)')
     .     'CHEMICAL_EROSION_BE_FAC_C=',
     .      chemical_erosion_be_fac_c,','
        endif
        write(hlp_frm,'(a,i2,a)') '(a,1p,',nstrat1,'(1x,g10.3,'',''))'
        do ib = 1, nstrat1
         if (crcstra1(ib).eq.'V' .and. volrecstart1(ib).ne.1.0e21_R8)
     >    write(34, '(a,i2,a,1pg10.3,a)')
     .    'VOLRECSTART(',ib,')=',volrecstart1(ib),','
        end do
        if (volrecinc.ne.0.0_R8) write(34, '(a,1pg10.3,a)')
     .   'VOLRECINC=',volrecinc,','
        write(34, '(a,1pg10.3,a)')
     .   'VOLRECWT=',volrecwt,','
        do ib = 1, 4
          if (rf_neut(ib).ne.1.0_R8) write(34,'(a,i1,a,1pg10.3,a)')
     .     'RF_NEUT(',ib,')=',rf_neut(ib)
        end do
        do ib = 1, nstrat1
          if (userfluxparm1(ib,1).ne.0.0_R8)
     >     write(34, '(a,i2,a,1pg10.3,a)')
     .     'USERFLUXPARM(',ib,',1)=',userfluxparm1(ib,1),','
          if (userfluxparm1(ib,2).ne.0.0_R8)
     >     write(34, '(a,i2,a,1pg10.3,a)')
     .     'USERFLUXPARM(',ib,',2)=',userfluxparm1(ib,2),','
        end do
        ib1 = 0
        do ib = 1, nstrai
          if (savestratum(ib).or.crcstra(ib).eq.'T') then
            ib1 = ib1+1
            do is = 1, nspecies1
              if (isp1(is).ne.0 .and. gpfc(is,ib).ne.0.0_R8)
     .         write(34,'(a,i3,a,i2,a,f5.2,a)')
     .          'GPFC(',isp1(is),',',ib1,')=',gpfc(is,ib),','
            end do
          end if
        end do
#ifndef NO_B2_CHEM_SPUT
        if (igass_chemical.ne.0 .and. igass_chemical.lt.ns)
     >   igass_chemical=iatm1(igass_chemical)
        if (igass_chemical.ne.0 .and. igass_chemical.lt.ns) then
          if (fchar_chemical.ne.0.0_R8) write(34, '(a,1pg10.3,a)')
     .     'FCHAR_CHEMICAL=',fchar_chemical,','
          if (igass_chemical.ne.0) write(34, '(a,i3,a)')
     .     'IGASS_CHEMICAL=',igass_chemical,','
          if (itsput_chemical.ne.0) write(34, '(a,i3,a)')
     .     'ITSPUT_CHEMICAL=',itsput_chemical,','
          if (issput_chemical.ne.0) write(34, '(a,i3,a)')
     .     'ISSPUT_CHEMICAL=',issput_chemical,','
          do ib = 0, nlim+nsts
            if (chemical_sputter_yield(ib).ne.0.0_R8)
     >       write(34,'(a,i3,a,1pg10.3,a)')
     .       'CHEMICAL_SPUTTER_YIELD(',ib,')=',
     .       chemical_sputter_yield(ib),','
          end do
        endif
#endif
        if (neutrals_time_mod.ne.0.0_R8) write(34,'(a,f10.7,a)')
     >   'NEUTRALS_TIME_MOD=', neutrals_time_mod,','
        if (neutrals_time_switch.ne.0.0_R8) write(34,'(a,f10.7,a)')
     >   'NEUTRALS_TIME_SWITCH=', neutrals_time_switch,','
        if (neutrals_time_mod.ne.0.0_R8) write(34,'(a,a,a)')
     >   'NEUTRALS_FILENAME=''', trim(neutrals_filename),''','
        write(hlp_frm,'(a,i2,a)') '(a,i1,a,',nstrat1,'(1x,i3,'',''))'
        if (maxval(lstrascl).gt.0) then
          do is = 1, nspecies1
            if (isp1(is).ne.0) then
              ib1 = 0
              do ib = 1, nstrat
                if (savestratum(ib).or.crcstra(ib).eq.'T') then
                  ib1 = ib1 + 1
                  lstrascl1(ib1,isp1(is))=lstrascl(ib,is)
                endif
              enddo
            endif
          enddo
          do is = 1, nspecies1
            if (isp1(is).ne.0) then
              write(34, hlp_frm) 'LSTRASCL(1,',isp1(is),')=',
     .         (lstrascl1(ib,isp1(is)),ib=1,nstrat1)
            end if
          end do
        endif
        if (l_neutrad.ne.0)
     >   write(34,'(a,i2,a)') 'L_NEUTRAD=',l_neutrad,','
        write(hlp_frm,'(a,i2,a)') '(a,',nmol,'(1x,i2,'',''))'
        if (maxval(lmolscl).gt.0)
     >   write(34,hlp_frm) 'LMOLSCL=',(lmolscl(is),is=1,nmol)
        if (maxval(lionscl).gt.0)
     >   write(34,hlp_frm) 'LIONSCL=',(lionscl(is),is=1,nion)
        write(hlp_frm,'(a,i3,a)') '(a,',ns1,'(1x,i3,'',''))'
        write(34,hlp_frm) 'B2ESPCR=',(b2espcr1(is),is=0,ns1-1)
        write(hlp_frm,'(a,i2,a)') '(a,',nspecies1,'(1x,i3,'',''))'
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
          write(hlp_frm,'(a,i2,a)') '(a,',ncns,'(1x,i3,'',''))'
          write(34,hlp_frm) 'LCNS=',(lcns(ib),ib=1,ncns)
        endif
        if (ntns.gt.0) then
          write(hlp_frm,'(a,i2,a)') '(a,',ntns,'(1x,i3,'',''))'
          write(34,hlp_frm) 'LTNS=',(ltns(ib),ib=1,ntns)
        endif
        nsns = 0
        do ib = 1, nstrai
          if (crcstra(ib).eq.'N'.or. crcstra(ib).eq.'W'.or.
     &        crcstra(ib).eq.'E') nsns = nsns + 1
          if (crcstra(ib).eq.'S') then
            if (rcpos1(ib).ne.-2) then
              if(mod(region1(rcstart1(ib),rcpos1(ib),0),4).ne.1) then
                nsns = nsns + 1
              end if
            else
              found = .true.
              do ix = 1, rc_list_size(ib)
                found = found .and.
     &           mod(region1(rc_list_x1(ix,ib),
     .                       rc_list_y1(ix,ib),0),4).ne.1
              end do
              if (found) nsns = nsns + 1
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
        write(34,'(a,1pg10.3,a)') 'NEUT_SCL_LIM=',neut_scl_lim,','
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
          write(hlp_frm,'(a,i2,a)')
     .     '(a,',n_spcsrf,'('''''''',a8,'''''''','',''))'
          write(34, hlp_frm) 'SPS_MTRL=',(sps_mtrl(ib),ib=1,n_spcsrf)
        endif
        do ib = 1, 100
          if (debug_flags(ib).ne.0) write(34,'(a,i3,a,i3,a)')
     .     'DEBUG_FLAGS(',ib,')=',debug_flags(ib),','
        end do
        if (dbg_eir_mc.ne.0) write(34,'(a,i3,a)')
     .   'DBG_EIR_MC=',dbg_eir_mc,','
        write(34,'(a)') '/'
        close(34)
        if (.not.streql(neutrals_filename,'b2.neutrals.parameters'))
     >   then
          new_filename =
     .     '../'//neutrals_filename(1:lnblnk(neutrals_filename))//'.2'
          inquire(file=trim(neutrals_filename),exist=todo)
          call read_neutrals_namelist(nx, ny, ns, .false.)
          call xertst (nstrai.le.nstraid,
     >      'Increase DEF_NSTRA in DIMENSIONS.F')
          inquire(file=trim(new_filename),exist=found)
          if (found) write(6,*) trim(new_filename), ' already present!'
        else
          todo = .false.
          found = .true.
        endif
        neutrals_namelist = 0
        use_eirene = 0
      enddo

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
*  ..initialize arrays used jointly
        call b2yt_wall_init(nncut1,nx1,ny1,nncut,nx,ny,nlim,nsts)
*  ..convert some wall arrays
        if(wall_namelist_used) then
          call b2yt_wall(nx1,ny1,nx,ny,ns,ns1,natm,natm1,nmol,nion,
     .     is_neutral1,iatm1,bundle_frac,zn1,l_chem_sput,
     .     itrack1,itrack2,itrack3)
*  ..write new wall namelist
          call b2yt_wall_final(nx1,ny1,ns1)
        else
*  ..in the absence of a namelist, put in defaults
          call get_wall_data(nx,ny,ns,BoRiS,wall_namelist_used)
        endif ! wall_namelist_used
*  ..read in the namelist data
        if(sputter_namelist_used) then
          open(99,file='../b2.sputter_save.parameters')
          read(99,sputter_save)
          close(99)
*  ..convert the sputter namelist data
          natmm = max(natm,natm1)
          nwallm = max(nwall,nwall1)
          allocate(sputter_tmp(nwall,0:nsm-1+natmm+nmol+nion))
          sputter_yield_yt = 0.0_R8
          sputter_yield2_yt = 0.0_R8
          sputter_tmp=0.0_R8
          do ix = -1, ny
            do iy = -1, ny
              if (xymap(ix,iy).gt.0) sputter_tmp(xymap(ix,iy),0:ns-1) =
     >         sputter_yield(ix,iy,0:ns-1,1)
            enddo
          enddo
          call b2yt_wall_bundle_sum(ns,ns1,nsm,is_neutral1,
     .     bundle_frac,natm,natm1,natmm,nmol,nion,iatm1,
     .     nwall,nwall1,nwallm,sputter_tmp,contr_area)
          do ib = 1, nbnd1
            sputter_yield_yt(imapx1(ib),imapy1(ib),0:ns1-1,1) =
     .       sputter_tmp(ib,0:ns1-1)
          end do
          sputter_tmp=0.0_R8
          do ix = -1, ny
            do iy = -1, ny
              if (xymap(ix,iy).gt.0) sputter_tmp(xymap(ix,iy),0:ns-1) =
     >         sputter_yield(ix,iy,0:ns-1,2)
            enddo
          enddo
          call b2yt_wall_bundle_sum(ns,ns1,nsm,is_neutral1,
     .     bundle_frac,natm,natm1,natmm,nmol,nion,iatm1,
     .     nwall,nwall1,nwallm,sputter_tmp,contr_area)
          do ib = 1, nbnd1
            sputter_yield_yt(imapx1(ib),imapy1(ib),0:ns1-1,2) =
     .       sputter_tmp(ib,0:ns1-1)
          end do
          sputter_tmp=0.0_R8
          do ix = -1, ny
            do iy = -1, ny
              if (xymap(ix,iy).gt.0) sputter_tmp(xymap(ix,iy),0:ns-1) =
     >         sputter_yield2(ix,iy,0:ns-1,1)
            enddo
          enddo
          call b2yt_wall_bundle_sum(ns,ns1,nsm,is_neutral1,
     .     bundle_frac,natm,natm1,natmm,nmol,nion,iatm1,
     .     nwall,nwall1,nwallm,sputter_tmp,contr_area)
          do ib = 1, nbnd1
            sputter_yield2_yt(imapx1(ib),imapy1(ib),0:ns1-1,1) =
     .       sputter_tmp(ib,0:ns1-1)
          end do
          sputter_tmp=0.0_R8
          do ix = -1, ny
            do iy = -1, ny
              if (xymap(ix,iy).gt.0) sputter_tmp(xymap(ix,iy),0:ns-1) =
     >         sputter_yield2(ix,iy,0:ns-1,2)
            enddo
          enddo
          call b2yt_wall_bundle_sum(ns,ns1,nsm,is_neutral1,
     .     bundle_frac,natm,natm1,natmm,nmol,nion,iatm1,
     .     nwall,nwall1,nwallm,sputter_tmp,contr_area)
          do ib = 1, nbnd1
            sputter_yield2_yt(imapx1(ib),imapy1(ib),0:ns1-1,2) =
     .       sputter_tmp(ib,0:ns1-1)
          end do
          deallocate(sputter_tmp)
*  ..write new sputtering namelist
          if(maxval(sputter_yield_yt).gt.0.0_R8.or.
     .       maxval(sputter_yield2_yt).gt.0.0_R8) then
            open(99,file='../b2.sputter_save.parameters.2')
            write(nout(5),'(a)') 'b2.sputter_save.parameters.2'
            write(6,*) 'Writing ../b2.sputter_save.parameters.2'
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
      endif

*  ..write new b2.feedback_save.parameters if needed
*   ..we can only do this after the species mapping has been done
      inquire(file='../b2.feedback_save.parameters',exist=todo)
      new_filename = '../b2.feedback_save.parameters.2'
      inquire(file=new_filename,exist=found)
      if (found) write(6,*) trim(new_filename),' already present!'
      if((lfeedback.or.b2stbc_feedback.ne.0.or.todo)
     &   .and. .not.found .and. .not.same_bundling) then
        call read_b2mod_feedback(99,'    ',feedback_namelist_used)
        if (maxval(saved_na_feedback_actuator).gt.0.0_R8) then
          todo = .false.
          is = 0
          do while (.not.todo .and. is.lt.nspecies)
            todo = todo .or.
     &        (saved_na_feedback_actuator(is).gt.0.0_R8 .and.
     &         ispecies1(is).ne.is)
            is = is+1
          end do
          if (todo) then
            allocate(saved_na_feedback_actuator1(0:natm-1))
            saved_na_feedback_actuator1 = 0.0_R8
            do is = 1, nspecies
              if (ispecies1(is).ne.-1)
     >         saved_na_feedback_actuator1(ispecies1(is)-1) =
     .          saved_na_feedback_actuator(is-1)
            end do
            write(*,*) 'Writing b2.feedback_save.parameters.2'
            write(nout(5),'(a)') trim(new_filename)
            open(99,file=new_filename)
            write(99,'(a9)') '&FEEDBACK'
            if (saved_cbshe_core.ne.0.0_R8)
     >       write(99,*) 'SAVED_CBSHE_CORE=',saved_cbshe_core,','
            if (saved_cbshi_core.ne.0.0_R8)
     >       write(99,*) 'SAVED_CBSHI_CORE=',saved_cbshi_core,','
            if (saved_cbsna_core.ne.0.0_R8)
     >       write(99,*) 'SAVED_CBSNA_CORE=',saved_cbsna_core,','
            if (saved_cbsch_core.ne.0.0_R8)
     >       write(99,*) 'SAVED_CBSCH_CORE=',saved_cbsch_core,','
            if (saved_cbsna_pfr1.ne.0.0_R8)
     >       write(99,*) 'SAVED_CBSNA_PFR1=',saved_cbsna_pfr1,','
            if (saved_cbsna_pfr2.ne.0.0_R8)
     >       write(99,*) 'SAVED_CBSNA_PFR2=',saved_cbsna_pfr2,','
            if (saved_cbsna_sol.ne.0.0_R8)
     >       write(99,*) 'SAVED_CBSNA_SOL= ',saved_cbsna_sol,','
            write(hlp_frm,'(a,i2,a)')
     >       '(a,1p,',nspecies1,'(1x,1e15.8,'',''))'
            write(99,hlp_frm) 'SAVED_NA_FEEDBACK_ACTUATOR=',
     >       saved_na_feedback_actuator1(0:nspecies1-1)
            write(99,'(a1)') '/'
            close(99)
            deallocate(saved_na_feedback_actuator1)
          end if
        end if
      end if

*   ..write new b2.feedback_control.parameters if needed
      inquire(file='../b2.feedback_control.parameters',exist=todo)
      new_filename='../b2.feedback_control.parameters.2'
      inquire(file=new_filename,exist=found)
      if (found) write(6,*) trim(new_filename),' already present!'
      if (todo.and..not.found.and.
     .   (.not.same_bundling.or.nx.ne.nx1.or.ny.ne.ny1)) then
        allocate(na_feedback_target1(0:natm-1))
        allocate(na_feedback_choice1(0:natm-1))
        allocate(na_feedback_option1(0:natm-1))
        allocate(na_feedback_actuator1(0:natm-1))
        allocate(na_feedback_ix11(0:natm-1))
        allocate(na_feedback_ix21(0:natm-1))
        allocate(na_feedback_iy11(0:natm-1))
        allocate(na_feedback_iy21(0:natm-1))
        allocate(na_feedback_ib1(0:natm-1))
        allocate(na_feedback_alpha1(0:natm-1))
        allocate(na_feedback_beta1(0:natm-1))
        allocate(na_feedback_const1(0:natm-1))
        allocate(na_feedback_time1(0:natm-1))
        allocate(na_feedback_puff_min1(0:natm-1))
        allocate(na_feedback_puff_max1(0:natm-1))
        allocate(na_feedback_overshoot1(0:natm-1))
        call read_b2mod_feedback_control(99)
        write(*,*) 'Writing b2.feedback_control.parameters.2'
        open(99,file=new_filename)
        write(nout(5),'(a)') trim(new_filename)
        write(99,'(a17)') '&FEEDBACK_CONTROL'
        if (vacuum_communication.ne.0) write(99,'(a,i2,a)')
     .   'VACUUM_COMMUNICATION=',vacuum_communication,','
        write (hlp_frm,'(a,i2,a)') '(a,',nvac,'(i2,'',''))'
        if (maxval(vacuum_communication_nreg).ne.0 .or.
     &      minval(vacuum_communication_nreg).ne.0) write(99,hlp_frm)
     .   'VACUUM_COMMUNICATION_NREG=',
     .   (vacuum_communication_nreg(ib),ib=1,nvac)
        if (maxval(vacuum_communication_method).ne.0 .or.
     &      minval(vacuum_communication_method).ne.0) write(99,hlp_frm)
     .   'VACUUM_COMMUNICATION_METHOD=',
     .   (vacuum_communication_method(ib),ib=1,nvac)
        write(hlp_frm,'(a,i3,a)')
     .   '(1x,a,i1,a,1p,',ns1,'(1x,g10.3,'',''))'
        do i2 = 1, nvac
          do i1 = 1, nvacreg
            iy = conformal_convert(non_commensurate,2,
     .       (vacuum_communication_ix1(i1,i2)+
     .        vacuum_communication_ix2(i1,i2))/2,
     .        vacuum_communication_iy(i1,i2),nx,ny,nx1,ny1,4)
            ix1 = conformal_convert(non_commensurate,1,
     .       vacuum_communication_ix1(i1,i2),
     .       vacuum_communication_iy(i1,i2),nx,ny,nx1,ny1,1)
            ix2 = conformal_convert(non_commensurate,1,
     .       vacuum_communication_ix2(i1,i2),
     .       vacuum_communication_iy(i1,i2),nx,ny,nx1,ny1,3)
            vacuum_communication_iy(i1,i2) = iy
            vacuum_communication_ix1(i1,i2) = ix1
            vacuum_communication_ix2(i1,i2) = ix2
            if (iy.ne.-2 .or. ix1.ne.-2 .or. ix2.ne.-2) then
              write(99,'(a,i1,a,i1,a,i4,a)')
     .         'VACUUM_COMMUNICATION_IY(',i1,',',i2,')=',
     .          vacuum_communication_iy(i1,i2),','
              write(99,'(a,i1,a,i1,a,i4,a)')
     .         'VACUUM_COMMUNICATION_IX1(',i1,',',i2,')=',
     .          vacuum_communication_ix1(i1,i2),','
              write(99,'(a,i1,a,i1,a,i4,a)')
     .         'VACUUM_COMMUNICATION_IX2(',i1,',',i2,')=',
     .          vacuum_communication_ix2(i1,i2),','
            end if
          end do
          if(maxval(vacuum_communication_alpha(0:ns-1,i2)).ne.0.0_R8.or.
     .       minval(vacuum_communication_alpha(0:ns-1,i2)).ne.0.0_R8)
     .     then
            cftmp(0:ns-1)=vacuum_communication_alpha(0:ns-1,i2)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            vacuum_communication_alpha(0:ns1-1,i2)=cftmp(0:ns1-1)
            write(99,hlp_frm) 'VACUUM_COMMUNICATION_ALPHA(0,',i2,')=',
     .       (vacuum_communication_alpha(is,i2),is=0,ns1-1)
          end if
          if(maxval(vacuum_communication_beta(0:ns-1,i2)).ne.1.0_R8.or.
     .       minval(vacuum_communication_beta(0:ns-1,i2)).ne.1.0_R8)
     .     then
            cftmp(0:ns-1)=vacuum_communication_beta(0:ns-1,i2)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            vacuum_communication_beta(0:ns1-1,i2)=cftmp(0:ns1-1)
            write(99,hlp_frm) 'VACUUM_COMMUNICATION_BETA(0,',i2,')=',
     .       (vacuum_communication_beta(is,i2),is=0,ns1-1)
          end if
        end do
        do is = 1, nspecies
          if (ispecies1(is).ne.-1) then
            na_feedback_target1(ispecies1(is)-1) =
     .       na_feedback_target(is-1)
            na_feedback_choice1(ispecies1(is)-1) =
     .       na_feedback_choice(is-1)
            na_feedback_option1(ispecies1(is)-1) =
     .       na_feedback_option(is-1)
            na_feedback_actuator1(ispecies1(is)-1) =
     .       na_feedback_actuator(is-1)
            na_feedback_ix11(ispecies1(is)-1) =
     .       conformal_convert(non_commensurate,1,
     .       na_feedback_ix1(is-1),na_feedback_iy1(is-1),
     .       nx,ny,nx1,ny1,1)
            na_feedback_iy11(ispecies1(is)-1) =
     .       conformal_convert(non_commensurate,2,
     .       na_feedback_ix1(is-1),na_feedback_iy1(is-1),
     .       nx,ny,nx1,ny1,0)
            na_feedback_ix21(ispecies1(is)-1) =
     .       conformal_convert(non_commensurate,1,
     .       na_feedback_ix2(is-1),na_feedback_iy2(is-1),
     .       nx,ny,nx1,ny1,3)
            na_feedback_iy21(ispecies1(is)-1) =
     .       conformal_convert(non_commensurate,2,
     .       na_feedback_ix2(is-1),na_feedback_iy2(is-1),
     .       nx,ny,nx1,ny1,2)
            na_feedback_ib1(ispecies1(is)-1) =
     .       na_feedback_ib(is-1)
            na_feedback_alpha1(ispecies1(is)-1) =
     .       na_feedback_alpha(is-1)
            na_feedback_beta1(ispecies1(is)-1) =
     .       na_feedback_beta(is-1)
            na_feedback_const1(ispecies1(is)-1) =
     .       na_feedback_const(is-1)
            na_feedback_puff_min1(ispecies1(is)-1) =
     .       na_feedback_puff_min(is-1)
            na_feedback_puff_max1(ispecies1(is)-1) =
     .       na_feedback_puff_max(is-1)
            na_feedback_overshoot1(ispecies1(is)-1) =
     .       na_feedback_overshoot(is-1)
            na_feedback_time1(ispecies1(is)-1) =
     .       na_feedback_time(is-1)
          end if
        end do
        write(hlp_frm,'(a,i2,a)') '(1x,a,',nspecies1,'(1x,i2,'',''))'
        if (maxval(na_feedback_choice1(0:nspecies1-1)).gt.0)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_CHOICE=  ',na_feedback_choice1(0:nspecies1-1)
        if (maxval(na_feedback_option1(0:nspecies1-1)).gt.0)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_OPTION=  ',na_feedback_option1(0:nspecies1-1)
        if (maxval(na_feedback_actuator1(0:nspecies1-1)).gt.0)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_ACTUATOR=',na_feedback_actuator1(0:nspecies1-1)
        write(hlp_frm,'(a,i2,a)') '(1x,a,',nspecies1,'(1x,i4,'',''))'
        if (maxval(na_feedback_ib1(0:nspecies1-1)).gt.0)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_IB= ',na_feedback_ib1(0:nspecies1-1)
        if (maxval(na_feedback_ix11(0:nspecies1-1)).gt.-2)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_IX1=',na_feedback_ix11(0:nspecies1-1)
        if (maxval(na_feedback_ix21(0:nspecies1-1)).gt.-2)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_IX2=',na_feedback_ix21(0:nspecies1-1)
        if (maxval(na_feedback_iy11(0:nspecies1-1)).gt.-2)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_IY1=',na_feedback_iy11(0:nspecies1-1)
        if (maxval(na_feedback_iy21(0:nspecies1-1)).gt.-2)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_IY2=',na_feedback_iy21(0:nspecies1-1)
        write(hlp_frm,'(a,i2,a)')
     .   '(1x,a,1p,',nspecies1,'(1x,g10.3,'',''))'
        if (maxval(na_feedback_target1(0:nspecies1-1)).gt.0.0_R8)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_TARGET=  ',na_feedback_target1(0:nspecies1-1)
        if (maxval(na_feedback_alpha1(0:nspecies1-1)).ne.0.001_R8 .or.
     &      minval(na_feedback_alpha1(0:nspecies1-1)).ne.0.001_R8)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_ALPHA=   ',na_feedback_alpha1(0:nspecies1-1)
        if (maxval(na_feedback_beta1(0:nspecies1-1)).ne.1.0_R8 .or.
     &      minval(na_feedback_beta1(0:nspecies1-1)).ne.1.0_R8)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_BETA=    ',na_feedback_beta1(0:nspecies1-1)
        if (maxval(na_feedback_const1(0:nspecies1-1)).ne.0.0_R8 .or.
     &      minval(na_feedback_const1(0:nspecies1-1)).ne.0.0_R8)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_CONST=   ',na_feedback_const1(0:nspecies1-1)
        if (maxval(na_feedback_puff_min1(0:nspecies1-1)).gt.0.0_R8)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_PUFF_MIN=',na_feedback_puff_min1(0:nspecies1-1)
        if (maxval(na_feedback_puff_max1(0:nspecies1-1)).gt.0.0_R8)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_PUFF_MAX=',na_feedback_puff_max1(0:nspecies1-1)
        if (maxval(na_feedback_overshoot1(0:nspecies1-1)).gt.0.0_R8)
     >   write(99,hlp_frm) 'NA_FEEDBACK_OVERSHOOT=',
     .     na_feedback_overshoot1(0:nspecies1-1)
        if (maxval(na_feedback_time1(0:nspecies1-1)).gt.0.0_R8)
     >   write(99,hlp_frm)
     .     'NA_FEEDBACK_TIME=',na_feedback_time1(0:nspecies1-1)
        write(99,'(a1)') '/'
        close(99)
        deallocate(na_feedback_target1)
        deallocate(na_feedback_choice1)
        deallocate(na_feedback_option1)
        deallocate(na_feedback_actuator1)
        deallocate(na_feedback_ix11)
        deallocate(na_feedback_ix21)
        deallocate(na_feedback_iy11)
        deallocate(na_feedback_iy21)
        deallocate(na_feedback_ib1)
        deallocate(na_feedback_alpha1)
        deallocate(na_feedback_beta1)
        deallocate(na_feedback_const1)
        deallocate(na_feedback_puff_min1)
        deallocate(na_feedback_puff_max1)
        deallocate(na_feedback_overshoot1)
        deallocate(na_feedback_time1)
      end if

*   ..we now do b2ah.dat, last because we need to overwrite zn, etc...
      if (.not.same_bundling.or.non_commensurate.eq.1) then
        write(*,*) 'Writing b2ah.dat.2'
        call cfopen (nout(3),'b2ah.dat.2','new','formatted')
        write(nout(5),'(a)') 'b2ah.dat.2'
        call b2ahx0 (ninp(5), idum(0))
        call xertst (idum(0).eq.ns,
     .   'incorrect species number in b2ah.dat')
        write(nout(3),'(a)') '*dimens         (ns; free format)'
        write(nout(3),*) ns1
        call b2rflb (ninp(5), lblcp(61:120))
        write(nout(3),'(a)')
     .   '*label          (lblcp, character*60; free format)'
        write(nout(3),'(a,a,a)') '''',trim(lblcp(61:120)),''''
        call b2rfcp (ninp(5), ns)
        zamin(0:ns1-1) = zamin1(0:ns1-1) ! need to do this transfer
        zamax(0:ns1-1) = zamax1(0:ns1-1) ! so that the write in
        zn(0:ns1-1) = zn1(0:ns1-1)       ! b2wfcp is correct
        am(0:ns1-1) = am1(0:ns1-1)
*    ..we now need to convert the cbrbrk values for non_commensurate cases
        if (non_commensurate.eq.1 .and. nncut.ge.1) then
          if (nncut.eq.1) then
            cbrbrk(cbirso) = rratio(leftcut1(1),nx1)
            cbrbrk(cbirso+1) = rratio(rightcut1(1),nx1)
          else if (nncut.eq.2) then
            cbrbrk(cbirso) = rratio(leftcut1(1),nx1)
            cbrbrk(cbirso+1) = rratio(rightcut1(1),nx1)
            cbrbrk(cbirso+2) = rratio(nxtl1,nx1)
            cbrbrk(cbirso+3) = rratio(rightcut1(2),nx1)
            cbrbrk(cbirso+4) = rratio(leftcut1(2),nx1)
            cbrbrk(cbirno+1) = rratio(nxtl1,nx1)
          end if
        end if
        do ireg = 0, cbnrso+cbnrno+cbnrwe+cbnrea-1
          do i = 0, 7
            cftmp(0:ns-1) = cbsna(i,0:ns-1,ireg)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            cbsna(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
            cftmp(0:ns-1) = cbsmo(i,0:ns-1,ireg)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            cbsmo(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
            cftmp(0:ns-1) = cbshi(i,0:ns-1,ireg)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            cbshi(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
            cftmp(0:ns-1) = cbrec(i,0:ns-1,ireg)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            cbrec(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
            cftmp(0:ns-1) = cbmsb(i,0:ns-1,ireg)
            call convert_bundle_r
     .       (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
            cbmsb(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
          enddo
        enddo
        do i = 0, 7
          cftmp(0:ns-1) = cfdf0(i,0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cfdf0(i,0:ns1-1) = cftmp(0:ns1-1)
          cftmp(0:ns-1) = cfdna(i,0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cfdna(i,0:ns1-1) = cftmp(0:ns1-1)
          cftmp(0:ns-1) = cfdpa(i,0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cfdpa(i,0:ns1-1) = cftmp(0:ns1-1)
          cftmp(0:ns-1) = cfvla(i,0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cfvla(i,0:ns1-1) = cftmp(0:ns1-1)
          cftmp(0:ns-1) = cfvsa(i,0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cfvsa(i,0:ns1-1) = cftmp(0:ns1-1)
          cftmp(0:ns-1) = cfhci(i,0:ns-1)
          call convert_bundle_r
     .     (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
          cfhci(i,0:ns1-1) = cftmp(0:ns1-1)
        enddo
        call b2wfcp (nout(3), 0, ns1)
        close(nout(3))
      endif

*   ..we now do b2ai.dat, using the data computed for b2ah.dat
      if (.not.same_bundling) then
        call cfopen (nout(4),'b2ai.dat.2','new','formatted')
        write(*,*) 'Writing b2ai.dat.2'
        write(nout(5),'(a)') 'b2ai.dat.2'
        call b2aix0 (ninp(7:6:-1), idum(0), idum(1), idum(2))
        call xertst (idum(2).eq.ns,
     .   'incorrect species number in b2ai.dat')
        write(nout(4),'(a)') '*dimens         (ns; free format)'
        write(nout(4),*) ns1
        call b2rflb (ninp(6), lblps(61:120))
        write(nout(4),'(a)')
     .   '*label          (lblps, character*60; free format)'
        write(nout(4),'(a,a,a)') '''',trim(lblps(61:120)),''''
*   ..read species data
        read (ninp(6),*) ligne
        do is = 0, ns-1
          read (ninp(6),*) id0, (zpdum(is,i),i=0,3)
        enddo
*   ..write species data
        write (nout(4),'(a)')
     &   '*specs     zamin     zamax     zn       am'
        do is = 0, ns1-1
          write (nout(4),fmt0) is,zamin(is),zamax(is),zn(is),am(is)
        enddo
*   ..read and write initial densities
        read (ninp(6),*) ligne
        write(nout(4),'(a)') trim(ligne)
        allocate(nbini(0:nsm-1))
        read (ninp(6),*) (nbini(is),is=0,ns-1)
        cftmp(0:ns-1) = nbini(0:ns-1)
        call convert_bundle_r(ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
        nbini(0:ns1-1) = cftmp(0:ns1-1)
        write (hlp_frm,'(a,i3,a)') '(1p,',ns1,'(1e9.2,1x))'
        write (nout(4),hlp_frm) (nbini(is),is=0,ns1-1)
*   ..read and write ttini
        read (ninp(6),*) ligne
        write(nout(4),'(a)') trim(ligne)
        read (ninp(6),*) ligne
        write(nout(4),'(a)') trim(ligne)
        close(nout(4))
        deallocate(nbini)
      endif
      deallocate(idum)

*   ..we now write a new SOLPS_FLUIDS file if needed
      if (.not.same_bundling) then
        write(*,*) 'Writing SOLPS_FLUIDS.2'
        allocate(species_name(0:ns1-1))
        do is = 0, ns1-1
          call species(is,species_name(is),.false.)
        enddo
        write(hlp_frm,'(a1,i3,a7)') '(',ns1,'(a,1x))'
        write(nout(5),'(a)') 'SOLPS_FLUIDS.2'
        open(unit=99,file='../SOLPS_FLUIDS.2')
        write(99,hlp_frm) (trim(species_name(is)),is=0,ns1-1)
        close(99)
        deallocate(species_name)
      endif

      call dealloc_b2mod_diag
      call dealloc_b2mod_wall
      call dealloc_b2mod_wall2
      call dealloc_b2mod_layer
      call dealloc_b2mod_ppout
      call dealloc_b2mod_b2plot
      call dealloc_b2mod_b2cmrc
      call dealloc_sputter_data
      call dealloc_b2mod_sputter
      call dealloc_input_profile
      call dealloc_b2mod_elements
      call dealloc_b2mod_boundary
      call dealloc_b2mod_neutrals
      call dealloc_b2mod_numerics
      call dealloc_b2mod_wall_init
      call dealloc_b2mod_wall2_init
      call dealloc_transport_namelist
      call dealloc_b2mod_transport_models
      deallocate(am1)
      deallocate(zn1)
      deallocate(zamin1)
      deallocate(zamax1)
      deallocate(is_neutral1)
      deallocate(na1)
      deallocate(ua1)
      deallocate(kinrgy1)
      deallocate(fna1)
      deallocate(uadia1)
      deallocate(sputter_yield_yt)
      deallocate(sputter_yield2_yt)
      deallocate(bundle_frac)
      deallocate(bundle_contr)
      deallocate(b2espcr1)
      deallocate(eb2ibcr)
      deallocate(ispecies1)
      deallocate(ifl1)
      deallocate(cftmp)
      deallocate(istmp)
      deallocate(iatm1)
      deallocate(ipls1)
#ifdef B25_EIRENE
      CALL EIRENE_DEALLOC_BRAEIR
      CALL EIRENE_DEALLOC_EIRBRA
      CALL EIRENE_DEALLOC_CCOUPL
      CALL EIRENE_EXTRAB25_CLEANUP
#endif
      
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
66666 FORMAT (I2,1X,A8,12(I3),1X,A10,1X,I2)

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

      if (i.eq.-2) then
        convert_position = -2
      else if (i.eq.-1 .and. iconf.le.1) then
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
      integer ixpos, iypos, iplg, nsize, nsize1
      integer eirleft(6), eirright(6)
      integer eirleft1(6), eirright1(6)
      integer nxtl, nxtr
      integer nxtl1, nxtr1
      logical atend
*   ..procedures
      integer conformal_convert, convert_position
      external conformal_convert, convert_position, get_nxt, b2yt_nxt,
     .         xerrab
cxpb Do not forget that in Eirene X is Radial and Y is poloidal !

      iypos = ix-1 ! Eirene IX is B2 IY
      eirleft = nx+7
      eirright= nx+7
      if (nncut.eq.0) then
        eirleft(1) = 1
        eirright(1)= nx+1
        eirleft1(1) = 1
        eirright1(1)= nx1+1
      else if (nncut.eq.1) then
        eirleft(1) = 1
        eirright(1)= leftcut(1)+1
        eirleft(2) = leftcut(1)+2
        eirright(2)= rightcut(1)+2
        eirleft(3) = rightcut(1)+3
        eirright(3)= nx+3
        eirleft1(1) = 1
        eirright1(1)= leftcut1(1)+1
        eirleft1(2) = leftcut1(1)+2
        eirright1(2)= rightcut1(1)+2
        eirleft1(3) = rightcut1(1)+3
        eirright1(3)= nx1+3
      else if (nncut.eq.2) then
        call get_nxt(nx, nxtl, nxtr)
        call b2yt_nxt(nx1, nxtl1, nxtr1)
        eirleft(1) = 1
        eirright(1)= leftcut(1)+1
        eirleft(2) = leftcut(1)+2
        eirright(2)= leftcut(2)+2
        eirleft(3) = leftcut(2)+3
        eirright(3)= nxtl+3
        eirleft(4) = nxtr+3+2
        eirright(4)= rightcut(2)+4
        eirleft(5) = rightcut(2)+5
        eirright(5)= rightcut(1)+5
        eirleft(6) = rightcut(1)+6
        eirright(6)= nx+6
        eirleft1(1) = 1
        eirright1(1)= leftcut1(1)+1
        eirleft1(2) = leftcut1(1)+2
        eirright1(2)= leftcut1(2)+2
        eirleft1(3) = leftcut1(2)+3
        eirright1(3)= nxtl1+3
        eirleft1(4) = nxtr1+3+2
        eirright1(4)= rightcut1(2)+4
        eirleft1(5) = rightcut1(2)+5
        eirright1(5)= rightcut1(1)+5
        eirleft1(6) = rightcut1(1)+6
        eirright1(6)= nx1+6
      endif
      iplg = 1
      do while (iy.gt.eirright(iplg).and.iplg.le.5)
        iplg = iplg+1
      enddo
      if (xory.eq.1) then
        if (ix.eq.0) then
          convert_eirene = 0
        else if (ix.gt.ny+1) then
          call xerrab ('Too large IX value in convert_eirene')
          convert_eirene = 0
        else
          ixpos = iy-1-iplg ! Eirene IY is B2 IX
          convert_eirene = 1 + conformal_convert
     .     (non_commensurate,2,ixpos,iypos,nx,ny,nx1,ny1,iside)
        end if
      else if (xory.eq.2) then
        if (iy.eq.0) then
          convert_eirene = 0
        elseif (iy.eq.eirleft(iplg)) then
          convert_eirene = eirleft1(iplg)
        elseif (iy.eq.eirright(iplg)) then
          convert_eirene = eirright1(iplg)
        elseif (iy.gt.nx+max(1,nncut*3)) then
          call xerrab ('Too large IY value in convert_eirene')
          convert_eirene = 0
        else
          ixpos = iy - eirleft(iplg)
          nsize = eirright(iplg) - eirleft(iplg)
          nsize1= eirright1(iplg)- eirleft1(iplg)
          atend = iside.eq.3
          convert_eirene = eirleft1(iplg) +
     .     convert_position(ixpos,nsize,nsize1,0,atend)
        endif
      else
        call xerrab ('In convert_eirene = 0 : xory test failed')
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
      integer ixiso, ixiso1, inseliy, inseliy1, nxtr, nxtl, nxtr1, nxtl1
      logical atend
*   ..procedures
      integer convert_position
      external convert_position, get_nxt, b2yt_nxt
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
       if (xory.eq.1 .and. nncut.ge.2) then ! treatment for DN cases X-direction
                                            ! depends on which side of the cut we are
         call get_nxt(nx, nxtl, nxtr)
         call b2yt_nxt(nx1, nxtl1, nxtr1)
         if (ix.le.nxtl) then
           nsize = nxtl
           nsize1 = nxtl1
           conformal_convert=convert_position(ipos,nsize,nsize1,0,atend)
         else if (ix.ge.nxtr) then
           ipos = ix - nxtr - 1
           nsize = nx - nxtr - 1
           nsize1 = nx1 - nxtr1 - 1
           conformal_convert= nxtr1 + 1 +
     .      convert_position(ipos,nsize,nsize1,0,atend)
         else
           conformal_convert = -2
         endif
       else
         conformal_convert=convert_position(ipos,nsize,nsize1,0,atend)
       endif
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
          conformal_convert=leftcut1(2)-1
         else if (ix.eq.leftcut(2)) then
          conformal_convert=leftcut1(2)
         else if (ix.eq.rightcut(1)-1) then
          conformal_convert=rightcut1(1)-1
         else if (ix.eq.rightcut(1)) then
          conformal_convert=rightcut1(1)
         else if (ix.eq.rightcut(2)-1) then
          conformal_convert=rightcut1(2)-1
         else if (ix.eq.rightcut(2)) then
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
     .            ix.gt.leftcut(1) .and. ix.lt.leftcut(2)-1) then
          nsize = leftcut(2)-leftcut(1)
          nsize1 = leftcut1(1)-leftcut1(1)
          conformal_convert = leftcut1(1) +
     .     convert_position(ix-leftcut(1),nsize,nsize1,3,atend)
         else if (nncut.eq.2 .and.
     .            ix.gt.leftcut(2) .and. ix.lt.ixiso-1) then
          nsize = ixiso-leftcut(2)
          nsize1 = ixiso1-leftcut1(2)
          conformal_convert = leftcut1(2) +
     .     convert_position(ix-leftcut(2),nsize,nsize1,2,atend)
         else if (nncut.eq.2 .and.
     .            ix.gt.ixiso+1 .and. ix.lt.rightcut(2)-1) then
          nsize = rightcut(2)-ixiso-1
          nsize1 = rightcut1(2)-ixiso1-1
          conformal_convert = ixiso1 + 1 +
     .     convert_position(ix-ixiso-1,nsize,nsize1,1,atend)
         else if (nncut.eq.2 .and.
     .            ix.gt.rightcut(2) .and. ix.lt.rightcut(1)-1) then
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
          write(*,*) 'ix conformal conversion not successful !'
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
          write(*,*) 'iy conformal conversion not successful !'
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
      integer is0(0:DEF_NSD-1)
      real (kind=R8) ::
     *   zn0(0:DEF_NSD-1), am0(0:DEF_NSD-1),
     *   zamin0(0:DEF_NSD-1), zamax0(0:DEF_NSD-1),
     *   bundle_contr(0:ns-1,0:DEF_NSD-1)
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
      integer is0(0:DEF_NSD-1)
      real (kind=R8) ::
     *   zn0(0:DEF_NSD-1), am0(0:DEF_NSD-1),
     *   zamin0(0:DEF_NSD-1), zamax0(0:DEF_NSD-1),
     *   bundle_contr(0:ns-1,0:DEF_NSD-1)
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
      logical found, unused(0:DEF_NSD-1)
*   ..procedures
      intrinsic nint
      external xerrab

cxpb We first need to find the total number of basic species nz
      unused = .true.
      ii = minloc(is0)
      nz = ii(1)-1
      do ic = 0, nz-1
        found = .false.
        iss = 0
        do while (.not.found .and. iss.lt.ns1)
          found = nint(zn1(iss)).eq.nint(zn0(ic)) .and.
     &     nint(zamin1(iss)).le.nint(zamin0(ic)) .and.
     &     nint(zamax1(iss)).ge.nint(zamax0(ic)) .and. unused(iss)
          if (.not.found) iss = iss + 1
        enddo
        if (iss.lt.ns1) then
          is1(is0(ic)) = iss
          am1(iss) = am0(ic)
          if (ic.lt.nz-1) then
            if (.not.(nint(zn1(iss)).eq.nint(zn0(ic+1)) .and.
     &          nint(zamin1(iss)).le.nint(zamin0(ic+1)) .and.
     &          nint(zamax1(iss)).ge.nint(zamax0(ic+1)))) then
              unused(iss) = .false.
            endif
          else
            unused(iss) = .false.
          endif
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

      subroutine average_bundle(ns1,ns,nsm,is_neutral1,bundle_frac,work)
      use b2mod_types
      use b2mod_b2cmpa
      implicit none
      integer ns1, ns, nsm
      integer is1, is, is0
      logical is_neutral1(0:ns1-1), same
      real (kind=R8) ::
     *   bundle_frac(0:ns1-1,0:ns-1), work(0:nsm-1),
     *   rold(0:ns-1), sum, value

cxpb We assign to the new bundle a weighted average of the values from the old bundle

      rold(0:ns-1)=work(0:ns-1)
      work(0:ns1-1)=0.0_R8
      do is1 = 0, ns1-1
        sum = 0.0_R8
        do is = 0, ns-1
          sum = sum + bundle_frac(is1,is)
          work(is1) = work(is1) + bundle_frac(is1,is)*rold(is)
        enddo
        if (sum.eq.0.0_R8) then  ! this is a new species
*  When all species have the same setting, we copy it unto the new species
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
        else  ! we do the average
          work(is1) = work(is1)/sum
        endif
      enddo

      return
      end subroutine average_bundle

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
*  When all species have the same setting, we copy it unto the new species
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

      subroutine e_bundle_transmit_r
     . (npls,npls1,e_bundle_contr,rold,rnew,itype)
      use b2mod_types
      implicit none
      integer npls, npls1, nplsm, itype
      integer ipls, ipls1
      real (kind=R8) :: e_bundle_contr(npls1,npls)
      real (kind=R8) :: rold(npls), rnew(npls1)
      real (kind=R8) :: sum
      logical :: done

cxpb We transmit the value of an NPLSI array to NPLS1
cxpb We must account for both the cases where bundles
cxpb are being expanded or are being consolidated
cxpb  itype = 1 : Simple sum
cxpb  itype = 2 : Weighted average
cxpb  itype = 3 : Inherit value from first old species

      rnew(1:npls1) = 0.0_R8
      do ipls1 = 1, npls1
        sum = 0.0_R8
        done = .false.
        do ipls = 1, npls
          if (e_bundle_contr(ipls1,ipls).ne.0.0_R8) then
            sum = sum + e_bundle_contr(ipls1,ipls)
            if (itype.ne.3) then
              rnew(ipls1) = rnew(ipls1) +
     .                      rold(ipls)*e_bundle_contr(ipls1,ipls)
            else if (itype.eq.3 .and. .not.done) then
              rnew(ipls1) = rold(ipls)
            endif
            done = .true.
          endif
        end do
        if (itype.eq.1) then  ! we simply sum the contributions
          continue
        else if (itype.eq.2) then  ! we average the contributions
          if (sum.ne.0.0_R8) then
            rnew(ipls1) = rnew(ipls1)/sum
          else
            rnew(ipls1) = 0.0_R8
          endif
        endif
      end do
      return
      end subroutine e_bundle_transmit_r

      subroutine e_bundle_transmit_i
     . (npls,npls1,e_bundle_contr,iold,inew,itype)
      use b2mod_types
      implicit none
      integer npls, npls1, nplsm, itype
      integer ipls, ipls1
      real (kind=R8) :: e_bundle_contr(npls1,npls)
      integer :: iold(npls), inew(npls1)
      integer :: sum
      logical :: done

cxpb We transmit the value of an NPLSI array to NPLS1
cxpb We must account for both the cases where bundles
cxpb are being expanded or are being consolidated
cxpb  itype = 1 : Simple sum
cxpb  itype = 2 : Average value
cxpb  itype = 3 : Inherit value from first old species

      inew(1:npls1) = 0
      do ipls1 = 1, npls1
        sum = 0
        done = .false.
        do ipls = 1, npls
          if (e_bundle_contr(ipls1,ipls).ne.0.0_R8) then
            sum = sum + 1
            if (itype.ne.3) then
              inew(ipls1) = inew(ipls1) + iold(ipls)
            else if (itype.eq.3 .and. .not.done) then
              inew(ipls1) = iold(ipls)
            endif
            done = .true.
          endif
        end do
        if (itype.eq.1) then  ! we simply sum the contributions
          continue
        else if (itype.eq.2) then  ! we average the contributions
          if (sum.ne.0.0_R8) then
            inew(ipls1) = inew(ipls1)/sum
          else
            inew(ipls1) = 0
          endif
        endif
      end do
      return
      end subroutine e_bundle_transmit_i

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

      integer function transform_tally(lvol,itally,ispc,iatm1,ipls1,
     .                                 nphoti,natmi,nmoli,nioni,nplsi,
     .                                 natmi1,nplsi1)
      implicit none
      logical lvol
      integer itally, ispc
      integer nphoti,natmi,natmi1,nmoli,nioni,nplsi,nplsi1
      integer iatm1(natmi), ipls1(nplsi)

      IF (ispc.eq.0) THEN
        transform_tally = 0
        return
      END IF

      IF (lvol) THEN    ! check volumetric tallies
        IF (itally.eq.1  .or. itally.eq.5  .or. itally.eq.10 .or.
     &      itally.eq.16 .or. itally.eq.22 .or. itally.eq.28 .or.
     &      itally.eq.63 .or. itally.eq.67 .or. itally.eq.71 .or.
     &      itally.eq.75 .or. itally.eq.85 .or. itally.eq.89 .or.
     &      itally.eq.93 .or. itally.eq.97) THEN ! these are NATM volume tallies
          transform_tally = iatm1(ispc)
        ELSE IF (itally.eq.14 .or. itally.eq.20 .or. itally.eq.26
     &      .or. itally.eq.32 .or. itally.eq.79 .or. itally.eq.-2
     &      .or. itally.eq.-4 .or. itally.eq.-5 .or. itally.eq.-6
     &      .or. itally.eq.-7 .or. itally.eq.-13) THEN ! these are NPLS volume tallies
          transform_tally = ipls1(ispc)
        ELSE
          transform_tally = ispc
        ENDIF
      ELSE              ! check surface tallies
        IF (itally.eq.1  .or. itally.eq.2  .or. itally.eq.3  .or.
     &      itally.eq.4  .or. itally.eq.5  .or. itally.eq.6  .or.
     &      itally.eq.26 .or. itally.eq.27 .or. itally.eq.28 .or.
     &      itally.eq.29 .or. itally.eq.30 .or. itally.eq.31 .or.
     &      itally.eq.51) THEN ! these are NATM surface tallies
          transform_tally = iatm1(ispc)
        ELSE IF (itally.eq.25 .or. itally.eq.50 .or. itally.eq.55) THEN ! these are NPLS surface tallies
          transform_tally = ipls1(ispc)
        ELSE IF (itally.eq.59) THEN ! these are NSPZ surface tallies
          IF (ispc.LE.NPHOTI) THEN
            transform_tally = ispc
          ELSE IF (ispc.LE.NPHOTI+NATMI) THEN
            transform_tally = NPHOTI+iatm1(ispc-NPHOTI)
          ELSE IF (ispc.LE.NPHOTI+NATMI+NMOLI+NIONI) THEN
            transform_tally = ispc-NATMI+NATMI1
          ELSE
            transform_tally = NPHOTI+NATMI1+NMOLI+NIONI+
     .             ipls1(ispc-(NPHOTI+NATMI+NMOLI+NIONI))
          END IF
        ELSE
          transform_tally = ispc
        ENDIF
      END IF

      return
      end function transform_tally
