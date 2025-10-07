      subroutine tallies(ns,tim,mpg,geo,state,state_ext,switch)
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_tallies
      use b2mod_switches
      use b2mod_constants
      use b2mod_subsys
      use b2us_map
      use b2us_geo
      use b2us_plasma
      implicit none
      integer, intent(in) :: ns
      real (kind=R8), intent(in) :: tim
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: geo
      type (B2State), intent (in) :: state
      type (B2StateExt), intent (in) :: state_ext
      type (switches), intent (in) :: switch
      real (kind=R8), save :: potave(0:nsdecl-1)
      real (kind=R8), allocatable :: tif(:), tnf(:), tef(:), rzf(:,:)
      real (kind=R8), allocatable :: ea(:,:), eaf(:,:), ptf(:,:)
#ifndef NO_CDF
#     include <netcdf.inc>
      integer status, ncid, xtypev, ndims, nvars, natts, unlimid, timd
      integer dimids(NF_MAX_VAR_DIMS), dimsize(NF_MAX_VAR_DIMS)
      integer, save :: ntstep
      real (kind=R8) :: tim_arr(1)
      real (kind=R8) :: naregtmp(0:mpg%nnreg(0),0:ns-1)
      real (kind=R8) :: cdftmp(0:mpg%nnreg(0))
      character*(NF_MAX_NAME) :: dimname(1), varname
      character*256 filename
      logical :: ex
      logical streql
      external streql
#endif
*   ..common blocks
      integer ncall
      integer iCv,iFc,is,ir,tallies_netcdf
      character*256 tally_format, species_format
      character*13, allocatable, save :: species_name(:)
      save ncall,tallies_netcdf
      data ncall/0/,tallies_netcdf/0/
#ifdef NO_LONG_LINES
      data tally_format, species_format
     *  /'(a,1p,8e12.4,/,(''+ '',10e12.4))',
     *   '(a,1x,i3,18(1x,a),/,(''+ '',20(1x,a)))'/
#else
      data tally_format, species_format
     *  /'(a,1p,(10000e12.4))','(a,1x,i3,999(1x,a))'/
#endif
* ..preliminaries
*   ..subprogram start-up calls
      call subini ('tallies')
*    ..initialise species_name
      if (ncall.eq.0) then
        call initialise_potentials(ns)
        allocate(species_name(0:ns-1))
        do is = 0, ns-1
          potave(is) = (potmin(is)+potmax(is))/2.0_R8
          call species(is,species_name(is),.false.)
        enddo
        call ipgeti('tallies_netcdf',tallies_netcdf)
      endif
*    ..produce tallies
      fnareg=0.0_R8
      fhpreg=0.0_R8
      fhmreg=0.0_R8
      fhireg=0.0_R8
      fhereg=0.0_R8
      fchreg=0.0_R8
      fhtreg=0.0_R8                                             !cjk
      fhjreg=0.0_R8
      nareg=0.0_R8
      tereg=0.0_R8
      nereg=0.0_R8
      ne2reg=0.0_R8
      tireg=0.0_R8
      nireg=0.0_R8
      poreg=0.0_R8
      volreg=0.0_R8
      qconvireg=0.0_R8                                          !cjk
      qconvereg=0.0_R8
      b2stbr_sna_reg=0.0_R8
      b2stbr_smo_reg=0.0_R8
      b2stbr_sne_reg=0.0_R8
      b2stbr_she_reg=0.0_R8
      b2stbr_shi_reg=0.0_R8
      b2stbr_sch_reg=0.0_R8
      b2stbc_sna_reg=0.0_R8
      b2stbc_smo_reg=0.0_R8
      b2stbc_sne_reg=0.0_R8
      b2stbc_she_reg=0.0_R8
      b2stbc_shi_reg=0.0_R8
      b2stbc_sch_reg=0.0_R8
      b2stbm_sna_reg=0.0_R8
      b2stbm_smo_reg=0.0_R8
      b2stbm_sne_reg=0.0_R8
      b2stbm_she_reg=0.0_R8
      b2stbm_shi_reg=0.0_R8
      b2stbm_sch_reg=0.0_R8
      b2sext_sna_reg=0.0_R8
      b2sext_smo_reg=0.0_R8
      b2sext_sne_reg=0.0_R8
      b2sext_she_reg=0.0_R8
      b2sext_shi_reg=0.0_R8
      b2sext_sch_reg=0.0_R8
      resco_reg=0.0_R8
      reshe_reg=0.0_R8
      reshi_reg=0.0_R8
      resmo_reg=0.0_R8
      resmt_reg=0.0_R8
      respo_reg=0.0_R8
      allocate(tif(1:mpg%nFc))
      allocate(tnf(1:mpg%nFc))
      allocate(tef(1:mpg%nFc))
      allocate(rzf(1:mpg%nFc,0:ns-1))
      call intface(mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,state%pl%ti,tif)
      call intface(mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,state%pl%tn,tnf)
      call intface(mpg%nCv,mpg%nFc,mpg%fcCv,geo%fcVol,state%pl%te,tef)
      call intface1
     .      (mpg%nCv,mpg%nFc,ns-1,mpg%fcCv,geo%fcVol,state%rt%rza,rzf)
      if (state_ext%ns.gt.0) then
        allocate(ea(1:mpg%nCv,0:state_ext%ns-1))
        allocate(eaf(1:mpg%nFc,0:state_ext%ns-1))
        allocate(ptf(1:mpg%nFc,0:state_ext%ns-1))
        do is = 0, state_ext%ns-1
          ea(:,is) = 0.5_R8*mp*state_ext%am(is)*state_ext%ua(:,is)**2+
     .               state_ext%ta(:,is)
        end do
        call intface1(mpg%nCv,mpg%nFc,ns-1,mpg%fcCv,geo%fcVol,ea,eaf)
        call intface1
     .     (mpg%nCv,mpg%nFc,ns-1,mpg%fcCv,geo%fcVol,state_ext%pt,ptf)
      end if
      do iFc = 1, mpg%nFc
        if (mpg%fcReg(iFc).eq.0) cycle
        do is = 0, ns-1
          fnareg(mpg%fcReg(iFc),is)=fnareg(mpg%fcReg(iFc),is)+
     1       state%dv%fna(iFc,0,is)+state%dv%fna(iFc,1,is)
          fhmreg(mpg%fcReg(iFc),is)=fhmreg(mpg%fcReg(iFc),is)+
     1       state%dv%fhm(iFc,0,is)+state%dv%fhm(iFc,1,is)
          fhpreg(mpg%fcReg(iFc),is)=fhpreg(mpg%fcReg(iFc),is)+
     1       state%dv%fhp(iFc,0,is)+state%dv%fhp(iFc,1,is)
          if (switch%tn_style.eq.0 .or. .not.is_neutral(is) .or.
     &       (switch%tn_style.eq.2 .and. nint(zn(is)).ne.1)) then
            qconvireg(mpg%fcReg(iFc),is)=                           !cjk
     1       qconvireg(mpg%fcReg(iFc),is)+tif(iFc)*
     1       (1.5_R8*(state%dv%fna_32(iFc,0,is)+
     1                state%dv%fna_32(iFc,1,is))+
     2        2.5_R8*(state%dv%fna_52(iFc,0,is)+
     2                state%dv%fna_52(iFc,1,is)))
          else if (switch%tn_style.eq.2 .and. nint(zn(is)).eq.1) then
            qconvireg(mpg%fcReg(iFc),is)=                           !cjk
     1       qconvireg(mpg%fcReg(iFc),is)+tnf(iFc)*
     1       (1.5_R8*(state%dv%fnn_32(iFc,0)+
     1                state%dv%fnn_32(iFc,1))+
     2        2.5_R8*(state%dv%fnn_52(iFc,0)+
     2                state%dv%fnn_52(iFc,1)))
          end if
          qconvereg(mpg%fcReg(iFc))=                                !cjk
     1     qconvereg(mpg%fcReg(iFc))+tef(iFc)*rzf(iFc,is)*
     1    (1.5_R8*(state%dv%fna_32(iFc,0,is)+state%dv%fna_32(iFc,1,is))+
     1     2.5_R8*(state%dv%fna_52(iFc,0,is)+state%dv%fna_52(iFc,1,is)))
        end do
        qconvereg(mpg%fcReg(iFc))=                                  !wuhsh
     1    qconvereg(mpg%fcReg(iFc))-tef(iFc)*
     1   (1.5_R8*state%dv%fch_p(iFc,0)/qe)
        fhireg(mpg%fcReg(iFc))=fhireg(mpg%fcReg(iFc))+
     1     state%dv%fhi(iFc,0)+state%dv%fhi(iFc,1)
        fhereg(mpg%fcReg(iFc))=fhereg(mpg%fcReg(iFc))+
     1     state%dv%fhe(iFc,0)+state%dv%fhe(iFc,1)
        fchreg(mpg%fcReg(iFc))=fchreg(mpg%fcReg(iFc))+
     1     state%dv%fch(iFc,0)+state%dv%fch(iFc,1)
        fhjreg(mpg%fcReg(iFc))=fhjreg(mpg%fcReg(iFc))+
     1     state%dv%fhj(iFc,0)+state%dv%fhj(iFc,1)
        fhtreg(mpg%fcReg(iFc))=fhtreg(mpg%fcReg(iFc))+
     1     state%dv%fht(iFc,0)+state%dv%fht(iFc,1)+
     2     state_ext%fhi(iFc,0)+state_ext%fhi(iFc,1)
        do is = 0, state_ext%ns-1
          fhtreg(mpg%fcReg(iFc))=fhtreg(mpg%fcReg(iFc))+        !cjk
     1     (1.0-switch%BoRiS)*eaf(iFc,is)*
     1     (state_ext%fa(iFc,0,is)+state_ext%fa(iFc,1,is))+
     2      ptf(iFc,is)*ev*
     2     (state_ext%fa(iFc,0,is)+state_ext%fa(iFc,1,is))
        enddo
      end do
      do iCv = 1, mpg%nCv
        if (mpg%cvReg(iCv).eq.0) cycle
        do is = 0, ns-1
          nareg(mpg%cvReg(iCv),is)=nareg(mpg%cvReg(iCv),is)+
     1       state%pl%na(iCv,is)*geo%cvVol(iCv)
          b2stbr_sna_reg(mpg%cvReg(iCv),is)=
     1       b2stbr_sna_reg(mpg%cvReg(iCv),is)+
     1       state%srw%b2stbr_sna(iCv,is)
          b2stbr_smo_reg(mpg%cvReg(iCv),is)=
     1       b2stbr_smo_reg(mpg%cvReg(iCv),is)+
     1       state%srw%b2stbr_smo(iCv,is)
          b2stbm_sna_reg(mpg%cvReg(iCv),is)=
     1       b2stbm_sna_reg(mpg%cvReg(iCv),is)+
     1       state%srw%b2stbm_sna(iCv,is)
          b2stbm_smo_reg(mpg%cvReg(iCv),is)=
     1       b2stbm_smo_reg(mpg%cvReg(iCv),is)+
     1       state%srw%b2stbm_smo(iCv,is)
          b2sext_sna_reg(mpg%cvReg(iCv),is)=
     1       b2sext_sna_reg(mpg%cvReg(iCv),is)+
     1       state_ext%sna(iCv,is)
          b2sext_smo_reg(mpg%cvReg(iCv),is)=
     1       b2sext_smo_reg(mpg%cvReg(iCv),is)+
     1       state_ext%smo(iCv,is)
          b2stbc_sna_reg(mpg%cvReg(iCv),is)=
     1       b2stbc_sna_reg(mpg%cvReg(iCv),is)+
     1       state%srw%b2stbc_sna(iCv,is)
          b2stbc_smo_reg(mpg%cvReg(iCv),is)=
     1       b2stbc_smo_reg(mpg%cvReg(iCv),is)+
     1       state%srw%b2stbc_smo(iCv,is)
          resco_reg(mpg%cvReg(iCv),is)=resco_reg(mpg%cvReg(iCv),is)+
     1       state%dv%resco(iCv,is)
          resmo_reg(mpg%cvReg(iCv),is)=resmo_reg(mpg%cvReg(iCv),is)+
     1       state%dv%resmo(iCv,is)
        end do
        tereg(mpg%cvReg(iCv))=tereg(mpg%cvReg(iCv))+
     1     state%pl%te(iCv)*state%dv%ne(iCv)*geo%cvVol(iCv)
        nereg(mpg%cvReg(iCv))=nereg(mpg%cvReg(iCv))+
     1     state%dv%ne(iCv)*geo%cvVol(iCv)
        ne2reg(mpg%cvReg(iCv))=ne2reg(mpg%cvReg(iCv))+
     1     state%dv%ne2(iCv)*geo%cvVol(iCv)
        if (switch%tn_style.eq.0) then
          tireg(mpg%cvReg(iCv))=tireg(mpg%cvReg(iCv))+
     1     state%pl%ti(iCv)*state%dv%ni(iCv,0)*geo%cvVol(iCv)
          nireg(mpg%cvReg(iCv))=nireg(mpg%cvReg(iCv))+
     1     state%dv%ni(iCv,0)*geo%cvVol(iCv)
        else if (switch%tn_style.eq.1) then
          tireg(mpg%cvReg(iCv))=tireg(mpg%cvReg(iCv))+
     1     state%pl%ti(iCv)*state%dv%ni(iCv,1)*geo%cvVol(iCv)
          nireg(mpg%cvReg(iCv))=nireg(mpg%cvReg(iCv))+
     1     state%dv%ni(iCv,1)*geo%cvVol(iCv)
        else if (switch%tn_style.eq.2) then
          tireg(mpg%cvReg(iCv))=tireg(mpg%cvReg(iCv))+
     1     state%pl%ti(iCv)*
     2    (state%dv%ni(iCv,0)-state%dv%nn(iCv))*geo%cvVol(iCv)
          nireg(mpg%cvReg(iCv))=nireg(mpg%cvReg(iCv))+
     1    (state%dv%ni(iCv,0)-state%dv%nn(iCv))*geo%cvVol(iCv)
        end if
        poreg(mpg%cvReg(iCv))=poreg(mpg%cvReg(iCv))+
     1     state%pl%po(iCv)*geo%cvVol(iCv)
        volreg(mpg%cvReg(iCv))=volreg(mpg%cvReg(iCv))+
     1     geo%cvVol(iCv)
        b2stbr_sne_reg(mpg%cvReg(iCv))=
     1     b2stbr_sne_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbr_sne(iCv)
        b2stbr_she_reg(mpg%cvReg(iCv))=
     1     b2stbr_she_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbr_she(iCv)
        b2stbr_shi_reg(mpg%cvReg(iCv))=
     1     b2stbr_shi_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbr_shi(iCv)
        b2stbr_sch_reg(mpg%cvReg(iCv))=
     1     b2stbr_sch_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbr_sch(iCv)
        b2stbc_sne_reg(mpg%cvReg(iCv))=
     1     b2stbc_sne_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbc_sne(iCv)
        b2stbc_she_reg(mpg%cvReg(iCv))=
     1     b2stbc_she_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbc_she(iCv)
        b2stbc_shi_reg(mpg%cvReg(iCv))=
     1     b2stbc_shi_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbc_shi(iCv)
        b2stbc_sch_reg(mpg%cvReg(iCv))=
     1     b2stbc_sch_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbc_sch(iCv)
        b2stbm_sne_reg(mpg%cvReg(iCv))=
     1     b2stbm_sne_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbm_sne(iCv)
        b2stbm_she_reg(mpg%cvReg(iCv))=
     1     b2stbm_she_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbm_she(iCv)
        b2stbm_shi_reg(mpg%cvReg(iCv))=
     1     b2stbm_shi_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbm_shi(iCv)
        b2stbm_sch_reg(mpg%cvReg(iCv))=
     1     b2stbm_sch_reg(mpg%cvReg(iCv))+
     1     state%srw%b2stbm_sch(iCv)
        b2sext_sne_reg(mpg%cvReg(iCv))=
     1     b2sext_sne_reg(mpg%cvReg(iCv))+
     1     state_ext%sne(iCv)
        b2sext_she_reg(mpg%cvReg(iCv))=
     1     b2sext_she_reg(mpg%cvReg(iCv))+
     1     state_ext%she(iCv)
        b2sext_shi_reg(mpg%cvReg(iCv))=
     1     b2sext_shi_reg(mpg%cvReg(iCv))+
     1     state_ext%shi(iCv)
        b2sext_sch_reg(mpg%cvReg(iCv))=
     1     b2sext_sch_reg(mpg%cvReg(iCv))+
     1     state_ext%sch(iCv)
        reshe_reg(mpg%cvReg(iCv))=reshe_reg(mpg%cvReg(iCv))+
     1     state%dv%reshe(iCv)
        reshi_reg(mpg%cvReg(iCv))=reshi_reg(mpg%cvReg(iCv))+
     1     state%dv%reshi(iCv)
        resmt_reg(mpg%cvReg(iCv))=resmt_reg(mpg%cvReg(iCv))+
     1     state%dv%resmt(iCv)
        respo_reg(mpg%cvReg(iCv))=respo_reg(mpg%cvReg(iCv))+
     1     state%dv%respo(iCv)
      end do
      do is = 0, ns-1
        rsanareg(0,is)=0.0_R8
        rsahireg(0,is)=0.0_R8
        rsamoreg(0,is)=0.0_R8
        rranareg(0,is)=0.0_R8
        rrahireg(0,is)=0.0_R8
        rramoreg(0,is)=0.0_R8
        rqahereg(0,is)=0.0_R8
        rqradreg(0,is)=0.0_R8
        rqbrmreg(0,is)=0.0_R8
        rcxnareg(0,is)=0.0_R8
        rcxhireg(0,is)=0.0_R8
        rcxmoreg(0,is)=0.0_R8
        b2stbr_sna_reg(0,is)=0.0_R8
        b2stbr_smo_reg(0,is)=0.0_R8
        b2stbc_sna_reg(0,is)=0.0_R8
        b2stbc_smo_reg(0,is)=0.0_R8
        b2stbm_sna_reg(0,is)=0.0_R8
        b2stbm_smo_reg(0,is)=0.0_R8
        b2sext_sna_reg(0,is)=0.0_R8
        b2sext_smo_reg(0,is)=0.0_R8
        nareg(0,is)=0.0_R8
        resco_reg(0,is)=0.0_R8
        resmo_reg(0,is)=0.0_R8
        do ir = 1, mpg%nnreg(0)
          rsanareg(0,is)=rsanareg(0,is)+rsanareg(ir,is)
          rsahireg(0,is)=rsahireg(0,is)+rsahireg(ir,is)
          rsamoreg(0,is)=rsamoreg(0,is)+rsamoreg(ir,is)
          rranareg(0,is)=rranareg(0,is)+rranareg(ir,is)
          rrahireg(0,is)=rrahireg(0,is)+rrahireg(ir,is)
          rramoreg(0,is)=rramoreg(0,is)+rramoreg(ir,is)
          rqahereg(0,is)=rqahereg(0,is)+rqahereg(ir,is)
          rqradreg(0,is)=rqradreg(0,is)+rqradreg(ir,is)
          rqbrmreg(0,is)=rqbrmreg(0,is)+rqbrmreg(ir,is)
          rcxnareg(0,is)=rcxnareg(0,is)+rcxnareg(ir,is)
          rcxhireg(0,is)=rcxhireg(0,is)+rcxhireg(ir,is)
          rcxmoreg(0,is)=rcxmoreg(0,is)+rcxmoreg(ir,is)
          nareg(0,is)=nareg(0,is)+nareg(ir,is)
          b2stbr_sna_reg(0,is)=b2stbr_sna_reg(0,is)+
     1     b2stbr_sna_reg(ir,is)
          b2stbr_smo_reg(0,is)=b2stbr_smo_reg(0,is)+
     1     b2stbr_smo_reg(ir,is)
          b2stbc_sna_reg(0,is)=b2stbc_sna_reg(0,is)+
     1     b2stbc_sna_reg(ir,is)
          b2stbc_smo_reg(0,is)=b2stbc_smo_reg(0,is)+
     1     b2stbc_smo_reg(ir,is)
          b2stbm_sna_reg(0,is)=b2stbm_sna_reg(0,is)+
     1     b2stbm_sna_reg(ir,is)
          b2stbm_smo_reg(0,is)=b2stbm_smo_reg(0,is)+
     1     b2stbm_smo_reg(ir,is)
          b2sext_sna_reg(0,is)=b2sext_sna_reg(0,is)+
     1     b2sext_sna_reg(ir,is)
          b2sext_smo_reg(0,is)=b2sext_smo_reg(0,is)+
     1     b2sext_smo_reg(ir,is)
          resco_reg(0,is)=resco_reg(0,is)+resco_reg(ir,is)
          resmo_reg(0,is)=resmo_reg(0,is)+resmo_reg(ir,is)
        enddo
      enddo
      do is = 0, ns-1
        fnareg(0,is)=0.0_R8
        fhpreg(0,is)=0.0_R8
        fhmreg(0,is)=0.0_R8
        qconvireg(0,is)=0.0_R8                                         !cjk
      enddo
      fhireg(0)=0.0_R8
      fhereg(0)=0.0_R8
      fhtreg(0)=0.0_R8
      fhjreg(0)=0.0_R8
      fchreg(0)=0.0_R8
      qconvereg(0)=0.0_R8                                              !cjk
      tereg(0)=0.0_R8
      nereg(0)=0.0_R8
      ne2reg(0)=0.0_R8
      tireg(0)=0.0_R8
      nireg(0)=0.0_R8
      poreg(0)=0.0_R8
      volreg(0)=0.0_R8
      b2stbr_sne_reg(0)=0.0_R8
      b2stbr_she_reg(0)=0.0_R8
      b2stbr_shi_reg(0)=0.0_R8
      b2stbr_sch_reg(0)=0.0_R8
      b2stbc_sne_reg(0)=0.0_R8
      b2stbc_she_reg(0)=0.0_R8
      b2stbc_shi_reg(0)=0.0_R8
      b2stbc_sch_reg(0)=0.0_R8
      b2stbm_sne_reg(0)=0.0_R8
      b2stbm_she_reg(0)=0.0_R8
      b2stbm_shi_reg(0)=0.0_R8
      b2stbm_sch_reg(0)=0.0_R8
      b2sext_sne_reg(0)=0.0_R8
      b2sext_she_reg(0)=0.0_R8
      b2sext_shi_reg(0)=0.0_R8
      b2sext_sch_reg(0)=0.0_R8
      b2divue(0)=0.0_R8
      b2divua(0)=0.0_R8
      b2exbe(0)=0.0_R8
      b2exba(0)=0.0_R8
      b2visa(0)=0.0_R8
      b2joule(0)=0.0_R8
      b2fraa(0)=0.0_R8
      b2str(0)=0.0_R8
      b2wrong1(0)=0.0_R8
      b2wrong2(0)=0.0_R8
      b2wrong3(0)=0.0_R8
      reshe_reg(0)=0.0_R8
      reshi_reg(0)=0.0_R8
      resmt_reg(0)=0.0_R8
      respo_reg(0)=0.0_R8
      rdneureg(0)=0.0_R8
      do ir = 1, mpg%nnreg(0)
        tereg(0)=tereg(0)+tereg(ir)
        nereg(0)=nereg(0)+nereg(ir)
        ne2reg(0)=ne2reg(0)+ne2reg(ir)
        tireg(0)=tireg(0)+tireg(ir)
        nireg(0)=nireg(0)+nireg(ir)
        poreg(0)=poreg(0)+poreg(ir)
        volreg(0)=volreg(0)+volreg(ir)
        b2stbr_sne_reg(0)=b2stbr_sne_reg(0)+b2stbr_sne_reg(ir)
        b2stbr_she_reg(0)=b2stbr_she_reg(0)+b2stbr_she_reg(ir)
        b2stbr_shi_reg(0)=b2stbr_shi_reg(0)+b2stbr_shi_reg(ir)
        b2stbr_sch_reg(0)=b2stbr_sch_reg(0)+b2stbr_sch_reg(ir)
        b2stbc_sne_reg(0)=b2stbc_sne_reg(0)+b2stbc_sne_reg(ir)
        b2stbc_she_reg(0)=b2stbc_she_reg(0)+b2stbc_she_reg(ir)
        b2stbc_shi_reg(0)=b2stbc_shi_reg(0)+b2stbc_shi_reg(ir)
        b2stbc_sch_reg(0)=b2stbc_sch_reg(0)+b2stbc_sch_reg(ir)
        b2stbm_sne_reg(0)=b2stbm_sne_reg(0)+b2stbm_sne_reg(ir)
        b2stbm_she_reg(0)=b2stbm_she_reg(0)+b2stbm_she_reg(ir)
        b2stbm_shi_reg(0)=b2stbm_shi_reg(0)+b2stbm_shi_reg(ir)
        b2stbm_sch_reg(0)=b2stbm_sch_reg(0)+b2stbm_sch_reg(ir)
        b2sext_sne_reg(0)=b2sext_sne_reg(0)+b2sext_sne_reg(ir)
        b2sext_she_reg(0)=b2sext_she_reg(0)+b2sext_she_reg(ir)
        b2sext_shi_reg(0)=b2sext_shi_reg(0)+b2sext_shi_reg(ir)
        b2sext_sch_reg(0)=b2sext_sch_reg(0)+b2sext_sch_reg(ir)
        b2qie(0)=b2qie(0)+b2qie(ir)
        b2divue(0)=b2divue(0)+b2divue(ir)
        b2divua(0)=b2divua(0)+b2divua(ir)
        b2exbe(0)=b2exbe(0)+b2exbe(ir)
        b2exba(0)=b2exba(0)+b2exba(ir)
        b2visa(0)=b2visa(0)+b2visa(ir)
        b2joule(0)=b2joule(0)+b2joule(ir)
        b2fraa(0)=b2fraa(0)+b2fraa(ir)
        b2str(0)=b2str(0)+b2str(ir)
        b2wrong1(0)=b2wrong1(0)+b2wrong1(ir)
        b2wrong2(0)=b2wrong2(0)+b2wrong2(ir)
        b2wrong3(0)=b2wrong3(0)+b2wrong3(ir)
        b2she(0)=b2she(0)+b2she(ir)
        b2shi(0)=b2shi(0)+b2shi(ir)
        b2she0(0)=b2she0(0)+b2she0(ir)
        b2shi0(0)=b2shi0(0)+b2shi0(ir)
        reshe_reg(0)=reshe_reg(0)+reshe_reg(ir)
        reshi_reg(0)=reshi_reg(0)+reshi_reg(ir)
        resmt_reg(0)=resmt_reg(0)+resmt_reg(ir)
        respo_reg(0)=respo_reg(0)+respo_reg(ir)
        rdneureg(0)=rdneureg(0)+rdneureg(ir)
      enddo
*    ..write tallies
      write(*,tally_format) 'Start tallies by region'
      write(*,'(a,2(1x,i2))') 'nnreg',mpg%nnreg
      write(*,species_format) 'species',ns,
     * (trim(species_name(is)),is=0,ns-1)
      if(ncall.eq.0) then
        write(*,tally_format) 'pot',potave(0:ns-1)
      endif
      write(*,tally_format) 'rsanareg',rsanareg
      write(*,tally_format) 'rsahireg',rsahireg
      write(*,tally_format) 'rsamoreg',rsamoreg
      write(*,tally_format) 'rranareg',rranareg
      write(*,tally_format) 'rrahireg',rrahireg
      write(*,tally_format) 'rramoreg',rramoreg
      write(*,tally_format) 'rqahereg',rqahereg
      write(*,tally_format) 'rqradreg',rqradreg
      write(*,tally_format) 'rqbrmreg',rqbrmreg
      write(*,tally_format) 'rcxnareg',rcxnareg
      write(*,tally_format) 'rcxhireg',rcxhireg
      write(*,tally_format) 'rcxmoreg',rcxmoreg
      write(*,tally_format) 'fnareg',fnareg
      write(*,tally_format) 'fhireg',fhireg
      write(*,tally_format) 'fhereg',fhereg
      write(*,tally_format) 'fhpreg',fhpreg
      write(*,tally_format) 'fhmreg',fhmreg
      write(*,tally_format) 'fchreg',fchreg
      write(*,tally_format) 'fhtreg',fhtreg                    !cjk
      write(*,tally_format) 'fhjreg',fhjreg
      write(*,tally_format) 'qconvireg',qconvireg
      write(*,tally_format) 'qconvereg',qconvereg
      write(*,tally_format) 'b2stbr_sna_reg',b2stbr_sna_reg
      write(*,tally_format) 'b2stbr_smo_reg',b2stbr_smo_reg
      write(*,tally_format) 'b2stbr_sne_reg',b2stbr_sne_reg
      write(*,tally_format) 'b2stbr_she_reg',b2stbr_she_reg
      write(*,tally_format) 'b2stbr_shi_reg',b2stbr_shi_reg
      write(*,tally_format) 'b2stbr_sch_reg',b2stbr_sch_reg
      write(*,tally_format) 'b2stbc_sna_reg',b2stbc_sna_reg
      write(*,tally_format) 'b2stbc_smo_reg',b2stbc_smo_reg
      write(*,tally_format) 'b2stbc_sne_reg',b2stbc_sne_reg
      write(*,tally_format) 'b2stbc_she_reg',b2stbc_she_reg
      write(*,tally_format) 'b2stbc_shi_reg',b2stbc_shi_reg
      write(*,tally_format) 'b2stbc_sch_reg',b2stbc_sch_reg
      write(*,tally_format) 'b2stbm_sna_reg',b2stbm_sna_reg
      write(*,tally_format) 'b2stbm_smo_reg',b2stbm_smo_reg
      write(*,tally_format) 'b2stbm_sne_reg',b2stbm_sne_reg
      write(*,tally_format) 'b2stbm_she_reg',b2stbm_she_reg
      write(*,tally_format) 'b2stbm_shi_reg',b2stbm_shi_reg
      write(*,tally_format) 'b2stbm_sch_reg',b2stbm_sch_reg
      write(*,tally_format) 'b2sext_sna_reg',b2sext_sna_reg
      write(*,tally_format) 'b2sext_smo_reg',b2sext_smo_reg
      write(*,tally_format) 'b2sext_sne_reg',b2sext_sne_reg
      write(*,tally_format) 'b2sext_she_reg',b2sext_she_reg
      write(*,tally_format) 'b2sext_shi_reg',b2sext_shi_reg
      write(*,tally_format) 'b2sext_sch_reg',b2sext_sch_reg
      write(*,tally_format) 'nareg',nareg
      write(*,tally_format) 'tereg',tereg
      write(*,tally_format) 'nereg',nereg
      write(*,tally_format) 'ne2reg',ne2reg
      write(*,tally_format) 'tireg',tireg
      write(*,tally_format) 'nireg',nireg
      write(*,tally_format) 'poreg',poreg
      write(*,tally_format) 'volreg',volreg
      write(*,tally_format) 'b2brem',b2brem
      write(*,tally_format) 'b2rad',b2rad
      write(*,tally_format) 'b2qie',b2qie
      write(*,tally_format) 'b2vdp',b2vdp
      write(*,tally_format) 'b2divue',b2divue
      write(*,tally_format) 'b2divua',b2divua
      write(*,tally_format) 'b2exbe',b2exbe
      write(*,tally_format) 'b2exba',b2exba
      write(*,tally_format) 'b2visa',b2visa
      write(*,tally_format) 'b2joule',b2joule
      write(*,tally_format) 'b2fraa',b2fraa
      write(*,tally_format) 'b2str',b2str
      write(*,tally_format) 'b2wrong1',b2wrong1
      write(*,tally_format) 'b2wrong2',b2wrong2
      write(*,tally_format) 'b2wrong3',b2wrong3
      write(*,tally_format) 'b2she',b2she
      write(*,tally_format) 'b2shi',b2shi
      write(*,tally_format) 'b2she0',b2she0
      write(*,tally_format) 'b2shi0',b2shi0
      write(*,tally_format) 'resco_reg',resco_reg
      write(*,tally_format) 'resmo_reg',resmo_reg
      write(*,tally_format) 'reshe_reg',reshe_reg
      write(*,tally_format) 'reshi_reg',reshi_reg
      write(*,tally_format) 'resmt_reg',resmt_reg
      write(*,tally_format) 'respo_reg',respo_reg
      write(*,tally_format) 'rdneureg',rdneureg
      write(*,tally_format) 'End tallies by region'
      deallocate(tef,tif,tnf,rzf)
      if (allocated(eaf)) deallocate(ea,eaf,ptf)

#ifndef NO_CDF
      if(tallies_netcdf.ne.0) then
        if(ncall.eq.0.or.mod(ncall,tallies_netcdf).eq.0) then
          write(*,tally_format) 'Start (CDF) tallies by region'
          filename='b2tallies.nc'
          call find_file(filename,ex)
          if(ncall.eq.0) then
            if (.not.ex.or.switch%b2mndr_stim.ge.0.0_R8) then
              ntstep = 0
              call define_cdf_b2mod_tallies(filename,ncid,mpg%nnreg,ns)
              status=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
              call check_cdf_status(status)
              write(6,'(a)') trim(filename)//' will be created'
            else if (ex.and.switch%b2mndr_stim.lt.0.0_R8) then
              status=nf_open(filename,NF_NOWRITE,ncid)
              call check_cdf_status(status)
              status=nf_inq(ncid,ndims,nvars,natts,unlimid)
              call check_cdf_status(status)
              status=nf_inq_varid(ncid,'times',timd)
              call check_cdf_status(status)
              status=nf_inq_var(ncid,timd,varname,xtypev,
     *                          ndims,dimids,natts)
              call check_cdf_status(status)
              call xertst(ndims.eq.1,
     *         'times field not one-dimensional in '//trim(filename))
              status=nf_inq_dim(ncid,dimids(1),dimname(1),dimsize(1))
              call check_cdf_status(status)
              ntstep = dimsize(1)
              status=nf_close(ncid)
              call check_cdf_status(status)
              status=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
              call check_cdf_status(status)
              write(6,'(a)') trim(filename)//' will be appended'
            else
              ntstep = 0
              write(6,'(a)') trim(filename)//' will be replaced'
              call define_cdf_b2mod_tallies(filename,ncid,mpg%nnreg,ns)
              status=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
              call check_cdf_status(status)
            end if
            write(6,'(a,i8)') trim(filename)//' : ntstep = ', ntstep
            if (ntstep.eq.0) call write_cdf_char_variable(ncid,-1,
     *       'species',species_name(0:ns-1))
          else
            status=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
            call check_cdf_status(status)
          endif
          ntstep=ntstep+1

          tim_arr=tim
          call write_cdf_variable(ncid,ntstep,
     *     'times',tim_arr)
          call write_cdf_variable(ncid,ntstep,
     *     'rsanareg',rsanareg)
          call write_cdf_variable(ncid,ntstep,
     *     'rsahireg',rsahireg)
          call write_cdf_variable(ncid,ntstep,
     *     'rsamoreg',rsamoreg)
          call write_cdf_variable(ncid,ntstep,
     *     'rranareg',rranareg)
          call write_cdf_variable(ncid,ntstep,
     *     'rrahireg',rrahireg)
          call write_cdf_variable(ncid,ntstep,
     *     'rramoreg',rramoreg)
          call write_cdf_variable(ncid,ntstep,
     *     'rqahereg',rqahereg)
          call write_cdf_variable(ncid,ntstep,
     *     'rqradreg',rqradreg)
          call write_cdf_variable(ncid,ntstep,
     *     'rqbrmreg',rqbrmreg)
          call write_cdf_variable(ncid,ntstep,
     *     'rcxnareg',rcxnareg)
          call write_cdf_variable(ncid,ntstep,
     *     'rcxhireg',rcxhireg)
          call write_cdf_variable(ncid,ntstep,
     *     'rcxmoreg',rcxmoreg)
          call write_cdf_variable(ncid,ntstep,
     *     'fnareg',fnareg)
          call write_cdf_variable(ncid,ntstep,
     *     'fhireg',fhireg)
          call write_cdf_variable(ncid,ntstep,
     *     'fhereg',fhereg)
          call write_cdf_variable(ncid,ntstep,
     *     'fhpreg',fhpreg)
          call write_cdf_variable(ncid,ntstep,
     *     'fhmreg',fhmreg)
          call write_cdf_variable(ncid,ntstep,
     *     'fchreg',fchreg)
          call write_cdf_variable(ncid,ntstep,                   !cjk
     *     'fhtreg',fhtreg)
          call write_cdf_variable(ncid,ntstep,
     *     'fhjreg',fhjreg)
          call write_cdf_variable(ncid,ntstep,
     *     'qconvireg',qconvireg)
          call write_cdf_variable(ncid,ntstep,
     *     'qconvereg',qconvereg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbr_sna_reg',b2stbr_sna_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbr_smo_reg',b2stbr_smo_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbr_sne_reg',b2stbr_sne_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbr_she_reg',b2stbr_she_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbr_shi_reg',b2stbr_shi_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbr_sch_reg',b2stbr_sch_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbc_sna_reg',b2stbc_sna_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbc_smo_reg',b2stbc_smo_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbc_sne_reg',b2stbc_sne_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbc_she_reg',b2stbc_she_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbc_shi_reg',b2stbc_shi_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbc_sch_reg',b2stbc_sch_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbm_sna_reg',b2stbm_sna_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbm_smo_reg',b2stbm_smo_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbm_sne_reg',b2stbm_sne_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbm_she_reg',b2stbm_she_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbm_shi_reg',b2stbm_shi_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2stbm_sch_reg',b2stbm_sch_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2sext_sna_reg',b2sext_sna_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2sext_smo_reg',b2sext_smo_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2sext_sne_reg',b2sext_sne_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2sext_she_reg',b2sext_she_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2sext_shi_reg',b2sext_shi_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2sext_sch_reg',b2sext_sch_reg)
          do is=0,ns-1
            naregtmp(:,is)=nareg(:,is)/volreg(:)
          enddo
          call write_cdf_variable(ncid,ntstep,
     *     'nareg',naregtmp)
          cdftmp=tereg/nereg/qe
          call write_cdf_variable(ncid,ntstep,
     *     'tereg',cdftmp)
          cdftmp=nereg/volreg
          call write_cdf_variable(ncid,ntstep,
     *     'nereg',cdftmp)
          cdftmp=ne2reg/volreg
          call write_cdf_variable(ncid,ntstep,
     *     'ne2reg',cdftmp)
          cdftmp=tireg/nireg/qe
          call write_cdf_variable(ncid,ntstep,
     *     'tireg',cdftmp)
          cdftmp=nireg/volreg
          call write_cdf_variable(ncid,ntstep,
     *     'nireg',cdftmp)
          cdftmp=poreg/volreg
          call write_cdf_variable(ncid,ntstep,
     *     'poreg',cdftmp)
          call write_cdf_variable(ncid,ntstep,
     *     'volreg',volreg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2bremreg',b2bremreg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2radreg',b2radreg)
          call write_cdf_variable(ncid,ntstep,
     *     'b2qie',b2qie)
!       call write_cdf_variable(ncid,ntstep,
!     *  'b2vdp',b2vdp)
          call write_cdf_variable(ncid,ntstep,
     *     'b2divue',b2divue)
          call write_cdf_variable(ncid,ntstep,
     *     'b2divua',b2divua)
          call write_cdf_variable(ncid,ntstep,
     *     'b2exbe',b2exbe)
          call write_cdf_variable(ncid,ntstep,
     *     'b2exba',b2exba)
          call write_cdf_variable(ncid,ntstep,
     *     'b2visa',b2visa)
          call write_cdf_variable(ncid,ntstep,
     *     'b2joule',b2joule)
          call write_cdf_variable(ncid,ntstep,
     *     'b2fraa',b2fraa)
          call write_cdf_variable(ncid,ntstep,
     *     'b2str',b2str)
          call write_cdf_variable(ncid,ntstep,
     *     'b2wrong1',b2wrong1)
          call write_cdf_variable(ncid,ntstep,
     *     'b2wrong2',b2wrong2)
          call write_cdf_variable(ncid,ntstep,
     *     'b2wrong3',b2wrong3)
          call write_cdf_variable(ncid,ntstep,
     *     'b2she',b2she)
          call write_cdf_variable(ncid,ntstep,
     *     'b2shi',b2shi)
          call write_cdf_variable(ncid,ntstep,
     *     'b2she0',b2she0)
          call write_cdf_variable(ncid,ntstep,
     *     'b2shi0',b2shi0)
          call write_cdf_variable(ncid,ntstep,
     *     'resco_reg',resco_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'resmo_reg',resmo_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'reshe_reg',reshe_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'reshi_reg',reshi_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'resmt_reg',resmt_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'respo_reg',respo_reg)
          call write_cdf_variable(ncid,ntstep,
     *     'rdneureg',rdneureg)
          write(*,tally_format) 'End (CDF) tallies by region'

          status=nf_close(ncid)
          call check_cdf_status(status)
        endif
      endif
#endif
* ..return
      ncall = ncall+1
      call subend ()
      return
      end subroutine tallies
