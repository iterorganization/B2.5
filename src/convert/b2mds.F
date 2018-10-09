      subroutine b2mds(nxd,nyd,nsd_local,in_jxi,in_jxa,in_jsep,BoRiS)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_anomalous_transport
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2plot
     * , only : textan, textmn, textin, textpl, textwl, textsp, textbl,
     * textcomponents, text0123, nx, ny, ns, jxa, jxi, jsep,
     * cr, cz, natmi, nmoli, nioni, wklng,
     * solpstop, alloc_b2mod_b2plot, dealloc_b2mod_b2plot
      use b2mod_wall
      use b2mod_work
      use b2mod_eirene_globals
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_transport_namelist
!      use b2mod_b2plot_debug
      use b2mod_elements
      use b2mod_sputter
!      use b2mod_movies
      use b2mod_external
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
#ifdef MDSPLUS
      use b2mod_mdsplus
#endif
      use b2mod_version
#ifdef B25_EIRENE
      use eirmod_cinit
     * , only : eirene_init_cinit
      use eirmod_clogau
      use eirmod_ctrcei
      use eirmod_comprt
      use eirmod_braeir
      use eirmod_eirbra
     * , only : eirene_alloc_eirbra, eirene_dealloc_eirbra
      use eirmod_extrab25
#else
      use b2mod_braeir
      use b2mod_eirbra
     & , smo_eir => smo
#endif
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_external
      implicit none
#     include <DIMENSIONS.F>
#ifdef MDSPLUS
#     include "mdslib.inc"
#endif
      integer nxd,nyd,nsd_local,in_jxi,in_jxa,in_jsep
      real (kind=R8) :: BoRiS
      character*24 hlp_frm
      integer iflag
      character*256 filename
#ifdef B25_EIRENE
      character*72 zeile
#endif
c
c  version : 27.07.96 20:37
c
      real (kind=R8) ::
     3   energy_balance(-1:nxd,-1:nyd),
     3   pr(-1:nxd,-1:nyd), energy_balance2(-1:nxd,-1:nyd),
     4   energy_balance_old(-1:nxd,-1:nyd)
      real (kind=R8) ::
     3   div_fna(-1:nxd,-1:nyd,0:nsd_local-1),
     3   div_fmo(-1:nxd,-1:nyd,0:nsd_local-1),
     3   div_fhm(-1:nxd,-1:nyd,0:nsd_local-1),
     3   div_fhp(-1:nxd,-1:nyd,0:nsd_local-1),
     3   div_fhet(-1:nxd,-1:nyd), div_fhit(-1:nxd,-1:nyd),
     4   div_fch(-1:nxd,-1:nyd), div_fhj(-1:nxd,-1:nyd)

      integer ispecies(0:nsd_local-1)

      integer is,iss,isp,nfla
c
      logical ex
c
      integer il,imi,ix,iy
#ifdef MDSPLUS
      character*512 textpl_cat, textwl_cat
      character*12 :: textdrift(0:1)
      integer :: status
      character :: mds_signal_name*(64), tdi*(172)
      real (kind=R8), allocatable ::
     *   vesl(:,:)
      integer i, icnt, nsw
      integer itmp(-1:nxd,-1:nyd)
      real (kind=R8) ::
     *   rdummy, rtmp, temp(-1:nxd,-1:nyd), hy1(-1:nxd,-1:nyd),
     5   tmp(-1:nxd,-1:nyd,0:nsd_local-1),
     5   cr_x(0:nxd,-1:nyd), cz_x(0:nxd,-1:nyd),
     6   cr_y(-1:nxd,0:nyd), cz_y(-1:nxd,0:nyd),
     7   flux_x(0:nxd,-1:nyd,0:nsd_local-1),
     7   flux_y(-1:nxd,0:nyd,0:nsd_local-1),
     8   dsrad(-1:nxd,-1:nyd),dspol(-1:nxd,-1:nyd),dspar(-1:nxd,-1:nyd),
     9   fhesep,fhisep,fnesep,fnisep,
     *   pot(0:nsd_local-1),poti(0:nsd_local-1),
     7   za(0:nsd_local-1)
#endif
c
      integer ismain, refl_model, sput_phys_model, use_eirene,
     . sput_phys_col, neutrals_namelist
      data ismain/1/, refl_model/1/, sput_phys_col/3/, use_eirene/0/,
     3     sput_phys_model/1/, neutrals_namelist/0/
      save ismain, refl_model, sput_phys_model,
     2     sput_phys_col, neutrals_namelist, use_eirene

#ifdef USE_PXFGETENV
      integer lenval, ierror
#else
#ifdef NAGFOR
      integer lenval, ierror
#endif
#endif
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
*     ------------------------------------------------------------------
      external subini, subend, ipgeti, xertst

c Function to check the correct parameters!
c
      call subini('b2mds')
#ifdef B25_EIRENE
      filename='input.dat'
      call find_file(filename,ex)
      if (ex) then
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
      if (.not.ex) then
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
        CALL EIRENE_SET_PARMMOD(1)
        CALL EIRENE_ALLOC_CLOGAU
        CALL EIRENE_INIT_CINIT
      endif
#endif
      NFLAI = NFL
#ifdef B25_EIRENE
      CALL EIRENE_ALLOC_BRAEIR(NDX,NDY,NFL,ifoff)
      CALL EIRENE_ALLOC_EIRBRA(NDX,NDY,NFL,NSTRA,ifoff)
      CALL EIRENE_ALLOC_CTRCEI
      CALL EIRENE_SETCON
      CALL EIRENE_EXTRAB25_WNEUINIT
#else
      call alloc_b2mod_braeir(nxd,nyd,nsd_local)
      call alloc_b2mod_eirbra(nxd,nyd,nsd_local,natm,npls,nstra)
#endif
      call alloc_b2mod_b2plot(nxd,nyd,nsd_local,ns_ext,wklng,
     1 natm,nmol,nion,nlimps)
#ifdef NO_GETENV
      solpstop=' '
#else
#ifdef NAGFOR
      call get_environment_variable('SOLPSTOP',
     .   status=ierror, length=lenval)
      if (ierror.eq.0) then
        call get_environment_variable('SOLPSTOP',value=solpstop)
      elseif (ierror.eq.1) then
        solpstop=' '
      endif
#else
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('SOLPSTOP', 0, solpstop, lenval, ierror)
#else
      call getenv ('SOLPSTOP', solpstop)
#endif
#endif
#endif
      cr=0.25_R8*sum(crx,3)
      cz=0.25_R8*sum(cry,3)
#ifdef MDSPLUS
      cr_x=0.0_R8
      cz_x=0.0_R8
      cr_y=0.0_R8
      cz_y=0.0_R8
      do iy=-1,ny
        do ix=0,nx
          cr_x(ix,iy)=(crx(ix,iy,0)+crx(ix,iy,2))/2.0_R8
          cz_x(ix,iy)=(cry(ix,iy,0)+cry(ix,iy,2))/2.0_R8
        enddo
      enddo
      do iy=0,ny
        do ix=-1,nx
          cr_y(ix,iy)=(crx(ix,iy,0)+crx(ix,iy,1))/2.0_R8
          cz_y(ix,iy)=(cry(ix,iy,0)+cry(ix,iy,1))/2.0_R8
        enddo
      enddo
#endif
c
      ns=max(nsd_local,1)
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',ns,'i3)'
      write(*,hlp_frm) 'zn: ',(nint(zn(is)),is=0,ns-1)
      do is=0,ns-1
        call species(is,textpl(is),.false.)
        call species(is,textwl(is),.false.)
      end do
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
      call ipgeti ('b2mndr_ismain', ismain)
      call xertst (0.le.ismain.and.ismain.lt.ns,
     &         'faulty b2mndr_ismain')
      call xertst (.not.is_neutral(ismain),
     &         'ismain cannot refer to a neutral species')
      textbl(0)=' '
      textcomponents(0)='x (poloidal)'
      textcomponents(1)='y'
      textcomponents(2)='z (toroidal)'
      textcomponents(3)='(total)'
      text0123(0)='0'
      text0123(1)='1'
      text0123(2)='2'
      text0123(3)='3'

cxpb  This call is moved up from 'f31' so that we use the original cell-faced
cxpb  quantities before their mapping to cell centers
      call b2xpni (nx, ny, ns, na, ni)
      call b2xppr (nx, ny, ns, ne, na, te, ti, pr)
      nfla=0
      do is=0,ns-1
        if(.not.is_neutral(is)) nfla=nfla+1
      enddo
      if(nfla.gt.NFL) then
        write(*,*) 'copy_background not called because nfla.gt.NFL',
     1   nfla,NFL
      elseif(nx.gt.NXDD.or.ny.gt.NYDD) then
        write(*,*)
     1   'copy_background not called because NX>NXDD or NY>NYDD',
     1   nx,nxdd,ny,nydd
      else
        call copy_background(nx,ny,ns,BoRiS)
      endif
c
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
        write(*,'(i3,1x,i3,1x,f5.1,1x,f5.1,1x,f5.1,4f11.3)')
     1   is,ispecies(is),zamin(is),zamax(is),zn(is),
     2   potimin(is),potimax(is),potmin(is),potmax(is)
      end do
      call alloc_b2mod_elements(ns, zn, am)
c calculate additional fluxes (fhj, fhm, fhp, fht)
      call b2xpen (nx, ny, ns,
     1 fna, fne, fch, fhe, fhi, kinrgy, rpt, te, ti, po, BoRiS,
     2 fhm, fhp, fhj, fht)
      if(allocated(fhi_ext)) fht=fht+fhi_ext
      do is=0,ns_ext-1
        do iy=-1,ny
          do ix=-1,nx
            if(leftix(ix,iy).ge.-1) then
              fht(ix,iy,0)=fht(ix,iy,0)+0.5_R8*fa_ext(ix,iy,0,is)*
     .          (1.0_R8-BoRiS)*(0.5_R8*am_ext(is)*mp*
     .          (ua_ext(leftix(ix,iy),leftiy(ix,iy),is)**2+
     .           ua_ext(ix,iy,is)**2)+
     .          (ta_ext(leftix(ix,iy),leftiy(ix,iy),is)+
     .           ta_ext(ix,iy,is)))
              fht(ix,iy,0)=fht(ix,iy,0)+0.5_R8*fa_ext(ix,iy,0,is)*
     .         (pt_ext(leftix(ix,iy),leftiy(ix,iy),is)+
     .          pt_ext(ix,iy,is))*ev
            endif
            if(bottomiy(ix,iy).ge.-1) then
              fht(ix,iy,1)=fht(ix,iy,1)+0.5_R8*fa_ext(ix,iy,1,is)*
     .         (1.0_R8-BoRiS)*(0.5_R8*am_ext(is)*mp*
     .         (ua_ext(bottomix(ix,iy),bottomiy(ix,iy),is)**2+
     .          ua_ext(ix,iy,is)**2)+
     .         (ta_ext(bottomix(ix,iy),bottomiy(ix,iy),is)+
     .          ta_ext(ix,iy,is)))
              fht(ix,iy,1)=fht(ix,iy,1)+0.5_R8*fa_ext(ix,iy,1,is)*
     .         (pt_ext(bottomix(ix,iy),bottomiy(ix,iy),is)+
     .          pt_ext(ix,iy,is))*ev
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
            div_fch(ix,iy)=
     1       +fch(rightix(ix,iy),rightiy(ix,iy),0)-fch(ix,iy,0)
     2       +fch(topix(ix,iy),topiy(ix,iy),1)-fch(ix,iy,1)
            div_fhj(ix,iy)=0
     1       +fhj(rightix(ix,iy),rightiy(ix,iy),0)-fhj(ix,iy,0)
     2       +fhj(topix(ix,iy),topiy(ix,iy),1)-fhj(ix,iy,1)
          else
            div_fhet(ix,iy)=0.0_R8
            div_fhit(ix,iy)=0.0_R8
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
c
cxpb  Bring in wall data and compute oft-needed wall quantities
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
      call ipgetr ('b2stbr_plate_temp', plate_temp)
      call xertst (0.0_R8.lt.plate_temp, 'faulty b2stbr_plate_temp')
      call ipgetr ('b2stbr_plate_thick', plate_thick)
      call xertst (0.0_R8.le.plate_thick, 'faulty b2stbr_plate_thick')
      call ipgeti ('b2stbr_sput_dst', sput_dst)
      call ipgeti ('b2stbr_sput_dst2', sput_dst2)
      call ipgeti ('b2stbr_sput_dst3', sput_dst3)
      call xertst (sput_dst.lt.ns, 'faulty b2stbr_sput_dst')
      call xertst (sput_dst2.lt.ns, 'faulty b2stbr_sput_dst2')
      call xertst (sput_dst3.lt.ns, 'faulty b2stbr_sput_dst3')
      call ipgeti ('b2stbr_sput_phys_col', sput_phys_col)
c the TRIM sputtering files have 9 columns corresponding to
c 0 15 30 45 55 65 75 80 85 degree impact
      call xertst (0.le.sput_phys_col.and.sput_phys_col.le.9,
     1 'faulty b2stbr_sput_phys_col')
      call ipgeti ('b2stbr_sput_phys_model', sput_phys_model)
      call alloc_b2mod_boundary(nx, ny, ns, nbcd, nnisod)
      call get_boundary_namelist
      if (boundary_namelist.eq.1)
     > call read_b2mod_boundary_namelist(nx, ny, ns, .true.)
      call xertst (nbc.le.nbcd, 'increase DEF_NBC in DIMENSIONS.F')
      call xertst (nniso.le.nnisod, 'increase DEF_NISO in DIMENSIONS.F')
      call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)
      call ipgeti ('b2mndr_eirene', use_eirene)
      call alloc_b2mod_neutrals
     > (nx, ny, ns, natm, nmol, nion, nstraid, nsts, sput_dst)
      if (neutrals_namelist.eq.1 .or. use_eirene.ne.0)
     > call read_neutrals_namelist(nx, ny, ns, .false.)
      call xertst (nstrat.le.nstraid,
     . 'Increase DEF_NSTRA in DIMENSIONS.F')
cxpb  Bring in wall data and compute oft-needed wall quantities
      write(*,*) 'Calling alloc_b2mod_sputter'
      call alloc_b2mod_sputter(nx,ny,ns)
      call alloc_b2mod_layer(nx,ny,nncut,nlim,nsts)
      call alloc_b2mod_wall_init(nx,ny,nncut)
      call alloc_b2mod_wall(ns,natm,nmol,nion)
      call init_wall(nx,ny)
      call read_b2mod_wall_namelist(99,'    ',ex)
      call xertst (ndepth_nml.gt.1.or.plate_model.eq.0,
     .  'Inconsistent values of ndepth_nml and plate_model')
      call get_wall_data(nx,ny,ns,BoRiS,ex)

#ifdef B25_EIRENE
#ifndef HPUX
      filename='fort.44'
#else
      filename='ftn44'
#endif
      call find_file(filename,ex)
#else
      ex = .false.
#endif
      if (.not.ex.or.use_eirene.eq.0) then
        natmi=0
        textan(0)=' ZERO'
        nmoli=0
        textmn(0)=' ZERO'
        nioni=0
        textin(0)=' ZERO'
        call eirene_mc_init(nx,ny)
      else
        open(unit=44,file=filename)
        call ngread(ismain,iflag)
        textwl(ns:ns+natmi-1)=textan(0:natmi-1)
        textwl(ns+natmi:ns+natmi+nmoli-1)=textmn(0:nmoli-1)
        textwl(ns+natmi+nmoli:ns+natmi+nmoli+nioni-1)=textin(0:nioni-1)
      endif

!XXX
      write(*,*) 'b2mds ...'
      if(in_jxi.ge.-1) jxi=in_jxi
      if(in_jxa.ge.-1) jxa=in_jxa
      if(in_jxi.ge.-1) jsep=in_jsep
    5 continue
      write(*,*) 'jxi, jxa, jsep : ', jxi, jxa, jsep

#ifdef MDSPLUS
      textdrift(0)='Diamagnetic'
      textdrift(1)='Radial'
      nsw=ns+natmi+nmoli+nioni

      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NX',nx+2)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NY',ny+2)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NS',ns)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NATM',natmi)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NMOL',nmoli)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NION',nioni)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NSW',nsw)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:NWALL',nbnd)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:TARG1',1)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:TARG2',-1)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:TARG3',-1)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:TARG4',nx)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:IMP',jxi+1)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:OMP',jxa+1)
      call mds_put_int('\TOP.SNAPSHOT.DIMENSIONS:SEP',jsep+1)

      tdi='build_with_units($,"m")'

      mds_signal_name='\TOP.SNAPSHOT.GRID:CR'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+2,0),cr,0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:CZ'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+2,0),cz,0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:CR_X'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx+1,ny+2,0),cr_x,0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:CZ_X'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2(idtype_double,nx+1,ny+2,0),cz_x,0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:CR_Y'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+1,0),cr_y,0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:CZ_Y'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+1,0),cz_y,0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:R'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+2,4,0),crx,
     * 0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      mds_signal_name='\TOP.SNAPSHOT.GRID:Z'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),
     * trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+2,4,0),cry,
     * 0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

      hy1(-1:nx,-1:ny)=hy(-1:nx,-1:ny)*qz(-1:nx,-1:ny,1)
      do ix=-1,nx
        dsrad(ix,jsep+1)=hy1(ix,jsep+1)/2.0_R8
        do iy=jsep+2,ny
          dsrad(ix,iy)=dsrad(ix,iy-1)+(hy1(ix,iy-1)+hy1(ix,iy))/2.0_R8
        enddo
        dsrad(ix,jsep)=-hy1(ix,jsep)/2.0_R8
        do iy=jsep-1,-1,-1
          dsrad(ix,iy)=dsrad(ix,iy+1)-(hy1(ix,iy+1)+hy1(ix,iy))/2.0_R8
        enddo
      enddo
      write(*,*) 'dsrad :',(dsrad(jxa,iy),iy=-1,ny)
      do iy=-1,ny
        dspol(-1,iy)=0.0_R8
        dspar(-1,iy)=0.0_R8
        do ix=0,nx
          dspol(ix,iy)=dspol(ix-1,iy)+(hx(ix-1,iy)+hx(ix,iy))/2.0_R8
          dspar(ix,iy)=dspar(ix-1,iy)+
     1     (hx(ix-1,iy)*abs(bb(ix-1,iy,2)/bb(ix-1,iy,0))+
     2      hx(ix,iy)*abs(bb(ix,iy,2)/bb(ix,iy,0)))/2.0_R8
        enddo
      enddo
      write(*,*) 'dspol :',(dspol(ix,jsep),ix=-1,nx)
      write(*,*) 'dspar :',(dspar(ix,jsep),ix=-1,nx)
      mds_signal_name='\TOP.SNAPSHOT.GRID:DSRAD'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+2,0),dsrad,0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)
      mds_signal_name='\TOP.SNAPSHOT.GRID:DSPOL'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+2,0),dspol,0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)
      mds_signal_name='\TOP.SNAPSHOT.GRID:DSPAR'
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nx+2,ny+2,0),dspar,0)
      write(*,*) 'mdsput2'
      call mds_err_msg(status)

c
      dbsolpsversion=b2fplasma_version
      call mds_put_string('\TOP.IDENT:SOLPSVERSION',dbsolpsversion)

c change dbversion below when MDS+ output format/data is amended/augmented
c 1: fht uses fhe and fhi instead of fhet and fhit
c 2: added fna{xy}_{35}2
c 3: fixed definition of \top.ident.impsep:kyi0 and \top.ident.ompsep:kyi0
c 4: added poti and pot
c 5: changed hcib and hci0 to include the radial neutral convective flux limits
c 6: added sxperp and qz
c 7: added new fields to identify b2f* file versions and svn version
c 8: changed from hy to hy1 in dsrad calculation
c 9: changed definition of za, pot and poti
c    added zamin, zamax, rza, rz2, rpt and rpi
c 10: added uadia, fch{xy}_{35}2, fch_p and kinrgy for recovery purposes
c 11: added fllimi & fllime
c 12: added ceqp
c 13: correctly save tmhacore, tmhasol & tmhadiv
c 14: includes contributions from external species in ne, pr, fht
c 15: correction of sign of electric energy flux
c 16: add \TOP.TIMEDEP.INPUT:DS \TOP.SNAPSHOT:FHJX \TOP.SNAPSHOT:FHJY
c     \TOP.TALLIES:fhjxreg \TOP.TALLIES:fhjyreg
c 17: corrected cr and cz locations for inner and outer separatrix midplanes
c 18: introduced new quantities from SOLPS_ITER
c     ENEUTRAD now with third dimension for atomic species
c 19: added GIT version
c 20: output na0, ua0, ne0, te0, ti0, po0
c 21: output smaf, smag, smav, smfr, smo, smpr, smpt, smq, smth, resco, reshe, reshi, resmo, resmt, respo
c 22: added various quantities for recovery of advanced SOLPS-ITER features
c     fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
c     fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
c     fni_32, fni_52, fne_32, fne_52,
c     fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
c     vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc
c 23: fht definition changed to 5/2 formulation
c 24: removed drift velocity contributions from kinrgy
c 25: debugged indexing error for eneutrad

      dbversion=25
      call mds_put_int('\TOP.IDENT:VERSION',dbversion)

      textpl_cat=textpl(0)
      do is=1,ns-1
        textpl_cat=trim(textpl_cat) // ' ' // textpl(is)
      enddo
      call mds_put_string_ns('\TOP.SNAPSHOT:TEXTPL',ns,textpl)
      call mds_put_string('\TOP.IDENT:SPECIES',textpl_cat)
      dbspecies=trim(textpl_cat)
c
      textwl_cat=textwl(0)
      do is=1,ns-1
        textwl_cat=trim(textwl_cat) // ' ' // textwl(is)
      enddo
      do is=ns,ns+natmi-1
        textwl_cat=trim(textwl_cat) // ' ' // textwl(is)
      enddo
      do is=ns+natmi,ns+natmi+nmoli-1
        textwl_cat=trim(textwl_cat) // ' ' // textwl(is)
      enddo
      do is=ns+natmi+nmoli,ns+natmi+nmoli+nioni-1
        textwl_cat=trim(textwl_cat) // ' ' // textwl(is)
      enddo
      call mds_put_string_ns('\TOP.SNAPSHOT:TEXTWL',
     1 ns+natmi+nmoli+nioni,textwl)
      call mds_put_string('\TOP.IDENT:WALL_SPECIES',textwl_cat)
c
      tdi='$'
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:ZAMIN','',
     1 ns,zamin(0:ns-1))
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:ZAMAX','',
     1 ns,zamax(0:ns-1))
      za(0:ns-1)=(zamin(0:ns-1)+zamax(0:ns-1))/2.0_R8
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:ZA','',ns,za(0:ns-1))
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:ZN','',ns,zn(0:ns-1))
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:AM','',ns,am(0:ns-1))
      poti(0:ns-1)=(potimin(0:ns-1)+potimax(0:ns-1))/2.0_R8
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:POTI','eV',ns,poti)
      pot(0:ns-1)=(potmin(0:ns-1)+potmax(0:ns-1))/2.0_R8
      call mds_put_double_ns('\TOP.SNAPSHOT.GRID:POT','eV',ns,pot)
      if(natmi.gt.0)
     1 call mds_put_string_ns('\TOP.SNAPSHOT:TEXTAN',natmi,textan)
      if(nmoli.gt.0)
     1 call mds_put_string_ns('\TOP.SNAPSHOT:TEXTMN',nmoli,textmn)
      if(nioni.gt.0)
     1 call mds_put_string_ns('\TOP.SNAPSHOT:TEXTIN',nioni,textin)
      call mds_put_string_ns('\TOP.SNAPSHOT:TEXTCOMP',4,textcomponents)
      call mds_put_string_ns('\TOP.SNAPSHOT:TEXTDRIFT',2,textdrift)
      call mds_put_string_ns('\TOP.SNAPSHOT:TEXTLIN',4,text0123)
c
      filename='mesh.extra'
      call find_file(filename,ex)
      if(ex) then
        open(99,file=filename,status='old')
        icnt=0
   10   continue
        read(99,*,end=20) (rdummy,ix=1,4)
        icnt=icnt+1
        goto 10
   20   continue
        write(*,*) 'Vessel size = ',icnt
        if (icnt.gt.0) then
          allocate(vesl(4,icnt))
          close(99)
          open(99,file='mesh.extra',status='old')
          do i=1,icnt
            read(99,*) (vesl(ix,i),ix=1,4)
          enddo
          close(99)
          tdi='build_with_units($,"m")'
          write(*,'(a,a,a,a)')
     1     'mdsput2 of ','\TOP.SNAPSHOT.GRID:VESSEL',' with ',trim(tdi)
          status=mdsput2(
     *     '\TOP.SNAPSHOT.GRID:VESSEL'//char(0),trim(tdi)//char(0),
     *     descr2(idtype_double,4,icnt,0),vesl,
     *     0)
          call mds_err_msg(status)
          deallocate(vesl)
        else
          write(*,*) 'Empty mesh.extra file found !'
        endif
      else
        write(*,*) 'No mesh.extra file found !'
      endif

      itmp=leftix+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:LEFTIX',nx+2,ny+2,itmp)
      itmp=leftiy+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:LEFTIY',nx+2,ny+2,itmp)
      itmp=rightix+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:RIGHTIX',nx+2,ny+2,itmp)
      itmp=rightiy+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:RIGHTIY',nx+2,ny+2,itmp)
      itmp=bottomix+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:BOTTOMIX',nx+2,ny+2,itmp)
      itmp=bottomiy+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:BOTTOMIY',nx+2,ny+2,itmp)
      itmp=topix+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:TOPIX',nx+2,ny+2,itmp)
      itmp=topiy+1
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:TOPIY',nx+2,ny+2,itmp)
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:REGVOL',nx+2,ny+2,region(:,:,0))
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:REGFLX',nx+2,ny+2,region(:,:,1))
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:REGFLY',nx+2,ny+2,region(:,:,2))
      call mds_put_int_rz(
     1 '\TOP.SNAPSHOT.GRID:XYMAP',nx+2,ny+2,xymap)
      call mds_put_int_ns(
     1 '\TOP.SNAPSHOT.GRID:IMAPX',nbnd,imapx)
      call mds_put_int_ns(
     1 '\TOP.SNAPSHOT.GRID:IMAPY',nbnd,imapy)

      temp=te/qe
      call mds_put_double_rz('\TOP.SNAPSHOT:TE','eV',nx+2,ny+2,temp,'')
      temp=te0/qe
      call mds_put_double_rz('\TOP.SNAPSHOT:TE0','eV',nx+2,ny+2,temp,'')
      temp=ti/qe
      call mds_put_double_rz('\TOP.SNAPSHOT:TI','eV',nx+2,ny+2,temp,'')
      temp=ti0/qe
      call mds_put_double_rz('\TOP.SNAPSHOT:TI0','eV',nx+2,ny+2,temp,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:NE','/m^3',nx+2,ny+2,ne,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:NE0','/m^3',
     &     nx+2,ny+2,ne0,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:PR','Pa',nx+2,ny+2,pr,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:PO','V',nx+2,ny+2,po,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:PO0','V',
     &     nx+2,ny+2,po0,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:VOL','m^3',nx+2,ny+2,vol,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:HX','m',nx+2,ny+2,hx,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:HY','m',nx+2,ny+2,hy,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:HY1','m',nx+2,ny+2,hy1,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:fllime','',
     1 nx+2,ny+2,fllime,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:fllimi','',
     1 nx+2,ny+2,fllimi,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:ceqp','',
     1 nx+2,ny+2,ceqp,'')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RZA','',
     1 nx+2,ny+2,ns,rza,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RZ2','',
     1 nx+2,ny+2,ns,rz2,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RPT','eV',
     1 nx+2,ny+2,ns,rpt,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RPI','eV',
     1 nx+2,ny+2,ns,rpi,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:NA','/m^3',
     1 nx+2,ny+2,ns,na,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:NA0','/m^3',
     1 nx+2,ny+2,ns,na0,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:UA','m/s',
     1 nx+2,ny+2,ns,ua,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:UA0','m/s',
     1 nx+2,ny+2,ns,ua0,'TEXTPL','')
      call mds_put_double_rzdns('\TOP.SNAPSHOT:UADIA','m/s',
     1 nx+2,ny+2,2,ns,uadia,'TEXTDRIFT','TEXTPL','')
      call mds_put_double_rzdns('\TOP.SNAPSHOT:VADIA','m/s',
     1 nx+2,ny+2,2,ns,vadia,'TEXTDRIFT','TEXTPL','')
      call mds_put_double_rzdns('\TOP.SNAPSHOT:WADIA','m/s',
     1 nx+2,ny+2,2,ns,wadia,'TEXTDRIFT','TEXTPL','')
      call mds_put_double_rzdns('\TOP.SNAPSHOT:VAECRB','m/s',
     1 nx+2,ny+2,2,ns,vaecrb,'TEXTDRIFT','TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:VEDIA','m/s',
     1 nx+2,ny+2,2,vedia,'TEXTDRIFT','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:VEECRB','m/s',
     1 nx+2,ny+2,2,veecrb,'TEXTDRIFT','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:KINRGY','J',
     1 nx+2,ny+2,ns,kinrgy,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:B','T',
     1 nx+2,ny+2,4,bb,'TEXTCOMP','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:QZ','',
     1 nx+2,ny+2,2,qz,'TEXTCOMP','')
      call mds_put_double_rz('\TOP.SNAPSHOT:BRHE','W',nx+2,ny+2,
     1 b2stbr_she,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:BRHI','W',nx+2,ny+2,
     1 b2stbr_shi,'')
      call mds_put_double_rzns('\TOP.SNAPSHOT:BRNA','particles/s',
     1 nx+2,ny+2,ns,b2stbr_sna,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:BRMO','N',
     1 nx+2,ny+2,ns,b2stbr_smo,'TEXTPL','')

      call mds_put_double_rzns('\TOP.SNAPSHOT:RQAHE','W',
     1 nx+2,ny+2,ns,rqahe,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RQRAD','W',
     1 nx+2,ny+2,ns,rqrad,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RQBRM','W',
     1 nx+2,ny+2,ns,rqbrm,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RSAHI','W',
     1 nx+2,ny+2,ns,rsahi,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RSANA','/s',
     1 nx+2,ny+2,ns,rsana,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RSAMO','m/s^2',
     1 nx+2,ny+2,ns,rsamo,'TEXTPL','')
       call mds_put_double_rzns('\TOP.SNAPSHOT:RRAHI','W',
     1 nx+2,ny+2,ns,rrahi,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RRANA','/s',
     1 nx+2,ny+2,ns,rrana,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RRAMO','m/s^2',
     1 nx+2,ny+2,ns,rramo,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RCXHI','W',
     1 nx+2,ny+2,ns,rcxhi,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RCXNA','/s',
     1 nx+2,ny+2,ns,rcxna,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RCXMO','m/s^2',
     1 nx+2,ny+2,ns,rcxmo,'TEXTPL','')

      call mds_put_double_rzns('\TOP.SNAPSHOT:D','m^2.s^-1',
     1 nx+2,ny+2,ns,dna0,'TEXTPL','')
      forall (ix=-1:nx,iy=-1:ny,is=0:ns-1)
     1 tmp(ix,iy,is)=dpa0(ix,iy,is)*(rza(ix,iy,is)*te(ix,iy)+ti(ix,iy))
      call mds_put_double_rzns('\TOP.SNAPSHOT:DP','m^2.s^-1',
     1 nx+2,ny+2,ns,tmp,'TEXTPL','')
      temp=hce0/ne
      call mds_put_double_rz('\TOP.SNAPSHOT:KYE','m^2.s^-1',
     1 nx+2,ny+2,temp,'')
      forall (ix=-1:nx,iy=-1:ny,is=0:ns-1)
     1 tmp(ix,iy,is)=hcib(ix,iy,is)/na(ix,iy,is)
      call mds_put_double_rzns('\TOP.SNAPSHOT:KYI','m^2.s^-1',
     1 nx+2,ny+2,ns,tmp,'TEXTPL','')
      temp=hci0(:,:)/ni(:,:,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:KYI0','m^2.s^-1',
     1 nx+2,ny+2,temp,'')
      tmp(:,:,:)=vla0(:,:,0,:)
      call mds_put_double_rzns('\TOP.SNAPSHOT:VLAX','m.s^-1',
     1 nx+2,ny+2,ns,tmp,'TEXTPL','')
      tmp(:,:,:)=vla0(:,:,1,:)
      call mds_put_double_rzns('\TOP.SNAPSHOT:VLAY','m.s^-1',
     1 nx+2,ny+2,ns,tmp,'TEXTPL','')
      forall (ix=-1:nx,iy=-1:ny,is=0:ns-1)
     1 tmp(ix,iy,is)=vsa0(ix,iy,is)/(mp*am(is)*na(ix,iy,is))
      call mds_put_double_rzns('\TOP.SNAPSHOT:VISC','m.s^-1',
     1 nx+2,ny+2,ns,tmp,'TEXTPL','')
      temp=sig0(:,:)/(qe*ne(:,:))
      call mds_put_double_rz('\TOP.SNAPSHOT:SIG','mho.m^-1',
     1 nx+2,ny+2,temp,'')
      temp=alf0(:,:)/(ne(:,:)*sqrt(qe/te(:,:)))
      call mds_put_double_rz('\TOP.SNAPSHOT:ALF','mho.m^-1',
     1 nx+2,ny+2,temp,'')

      flux_x(0:nx,-1:ny,0:ns-1)=fna(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAX','/s',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fna(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAY','/s',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fna_32(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAX_32','/s',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fna_32(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAY_32','/s',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fna_52(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAX_52','/s',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fna_52(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAY_52','/s',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fna_mdf(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAX_MDF','/s',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fna_mdf(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAY_MDF','/s',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fna_fcor(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAX_FCOR','/s',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fna_fcor(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAY_FCOR','/s',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fna_nodrift(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAX_NODRIFT','/s',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fna_nodrift(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAY_NODRIFT','/s',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fna_he(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAX_HE','/s',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fna_he(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAY_HE','/s',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fnaPSch(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAPSCHX','/s',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fnaPSch(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAPSCHY','/s',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fna_eir(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAX_EIR','/s',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fna_eir(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FNAY_EIR','/s',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fmo(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FMOX','N',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fmo(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FMOY','N',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fhm(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FHMX','W',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fhm(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FHMY','W',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0:ns-1)=fhp(0:nx,-1:ny,0,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FHPX','W',
     1 nx+1,ny+2,ns,flux_x,'TEXTPL','_X')
      flux_y(-1:nx,0:ny,0:ns-1)=fhp(-1:nx,0:ny,1,0:ns-1)
      call mds_put_double_rzns('\TOP.SNAPSHOT:FHPY','W',
     1 nx+2,ny+1,ns,flux_y,'TEXTPL','_Y')
      flux_x(0:nx,-1:ny,0)=fne_eir(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FNEX_EIR','/s',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fne_eir(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FNEY_EIR','/s',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fne_32(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FNEX_32','/s',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fne_32(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FNEY_32','/s',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fne_52(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FNEX_52','/s',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fne_52(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FNEY_52','/s',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fni_32(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FNIX_32','/s',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fni_32(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FNIY_32','/s',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fni_52(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FNIX_52','/s',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fni_52(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FNIY_52','/s',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fhe(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHEX','W',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fhe(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHEY','W',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fhe_mdf(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHEX_MDF','W',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fhe_mdf(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHEY_MDF','W',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fhePSch(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHEPSCHX','W',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fhePSch(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHEPSCHY','W',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fhe_eir(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHEX_EIR','W',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fhe_eir(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHEY_EIR','W',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fhi(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHIX','W',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fhi(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHIY','W',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fhi_mdf(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHIX_MDF','W',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fhi_mdf(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHIY_MDF','W',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fhiPSch(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHIPSCHX','W',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fhiPSch(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHIPSCHY','W',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fhi_eir(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHIX_EIR','W',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fhi_eir(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHIY_EIR','W',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fch(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHX','A',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fch(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHY','A',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fch_32(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHX_32','A',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fch_32(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHY_32','A',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fch_52(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHX_52','A',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fch_52(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHY_52','A',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fch_p(0:nx,-1:ny,0)             !srv 10.06.09
      call mds_put_double_rz('\TOP.SNAPSHOT:FCH_P','A',
     1 nx+1,ny+2,flux_x,'_X')
      flux_x(0:nx,-1:ny,0)=fchdia(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHDIAX','A',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fchdia(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHDIAY','A',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fchin(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHINX','A',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fchin(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHINY','A',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fchvispar(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHVISPARX','A',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fchvispar(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHVISPARY','A',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fchvisper(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHVISPERX','A',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fchvisper(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHVISPERY','A',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fchvisq(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHVISQX','A',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fchvisq(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHVISQY','A',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fchinert(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHINERTX','A',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fchinert(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FCHINERTY','A',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fht(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHTX','W',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fht(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHTY','W',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=fhj(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHJX','W',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=fhj(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FHJY','W',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=floe_noc(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FLOEX_NOC','/s',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=floe_noc(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FLOEY_NOC','/s',
     1 nx+2,ny+1,flux_y,'_Y')
      flux_x(0:nx,-1:ny,0)=floi_noc(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:FLOIX_NOC','/s',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=floi_noc(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:FLOIY_NOC','/s',
     1 nx+2,ny+1,flux_y,'_Y')

      flux_x(0:nx,-1:ny,0)=gs(0:nx,-1:ny,0)
      call mds_put_double_rz('\TOP.SNAPSHOT:SX','m^2',
     1 nx+1,ny+2,flux_x,'_X')
      do iy=-1,ny
        do ix=0,nx
          flux_x(ix,iy,0)=gs(ix,iy,0)*
     1     (vol(ix,iy)*qz(leftix(ix,iy),leftiy(ix,iy),1)+
     2     vol(leftix(ix,iy),leftiy(ix,iy))*qz(ix,iy,1))/
     3     (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
        enddo
      enddo
      call mds_put_double_rz('\TOP.SNAPSHOT:SXPERP','m^2',
     1 nx+1,ny+2,flux_x,'_X')
      flux_y(-1:nx,0:ny,0)=gs(-1:nx,0:ny,1)
      call mds_put_double_rz('\TOP.SNAPSHOT:SY','m^2',
     1 nx+2,ny+1,flux_y,'_Y')

      tdi='build_with_units($,"K")'
      mds_signal_name='\TOP.SNAPSHOT:WALL_TEMP'
      status=mdsput2(
     * trim(mds_signal_name)//char(0),trim(tdi)//char(0),
     * descr2(idtype_double,nbnd,ndepth_nml,0),target_temp,0)
      write(*,*) 'mdsput2 of ',trim(mds_signal_name),' with ',trim(tdi)
      call mds_err_msg(status)

      call mds_put_double_rzdns('\TOP.SNAPSHOT:SMAF','N',
     1 nx+2,ny+2,4,ns,b2npmo_smaf,'TEXTLIN','TEXTPL','')
      call mds_put_double_rzdns('\TOP.SNAPSHOT:SMAG','N',
     1 nx+2,ny+2,4,ns,b2npmo_smag,'TEXTLIN','TEXTPL','')
      call mds_put_double_rzdns('\TOP.SNAPSHOT:SMAV','N',
     1 nx+2,ny+2,4,ns,b2npmo_smav,'TEXTLIN','TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:SMFR','N',
     1 nx+2,ny+1,ns,smfr,'TEXTPL','')
      call mds_put_double_rzdns('\TOP.SNAPSHOT:SMO','N',
     1 nx+2,ny+2,4,ns,smo,'TEXTLIN','TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:SMPR','N',
     1 nx+2,ny+1,ns,smpr,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:SMPT','N',
     1 nx+2,ny+1,ns,smpt,'TEXTPL','')
      call mds_put_double_rzdns('\TOP.SNAPSHOT:SMQ','N',
     1 nx+2,ny+2,4,ns,smq,'TEXTLIN','TEXTPL','')
      tmp=smpr+smpt
      call mds_put_double_rzns('\TOP.SNAPSHOT:SMTH','N',
     1 nx+2,ny+1,ns,tmp,'TEXTPL','')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RESCO','#/s',
     1 nx+2,ny+1,ns,resco,'TEXTPL','')
      call mds_put_double_rz('\TOP.SNAPSHOT:RESHE','W',
     1 nx+2,ny+2,reshe,'')
      call mds_put_double_rz('\TOP.SNAPSHOT:RESHI','W',
     1 nx+2,ny+2,reshi,'')
      call mds_put_double_rzns('\TOP.SNAPSHOT:RESMO','N',
     1 nx+2,ny+1,ns,resmo,'TEXTPL','')
      call mds_put_double_rz('\TOP.SNAPSHOT:RESMT','N',
     1 nx+2,ny+2,resmt,'')
       call mds_put_double_rz('\TOP.SNAPSHOT:RESPO','V',
     1 nx+2,ny+2,respo,'')

#ifdef B25_EIRENE
      if(natmi.gt.0) then
        forall (ix=-1:nx,iy=-1:ny,is=0:natmi-1)
     1   tmp(ix,iy,is)=dab2(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:DAB2','/m^3',nx+2,
     1   ny+2,natmi,tmp,'TEXTAN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:natmi-1)
     1   tmp(ix,iy,is)=tab2(ix+1,iy+1,is+1,1)/ev
        call mds_put_double_rzns('\TOP.SNAPSHOT:TAB2','eV',nx+2,
     1   ny+2,natmi,tmp,'TEXTAN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:natmi-1)
     1   tmp(ix,iy,is)=rfluxa(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:RFLA','s^-1.m^-2',nx+2,
     1   ny+2,natmi,tmp,'TEXTAN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:natmi-1)
     1   tmp(ix,iy,is)=refluxa(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:REFA','W.m^-2',nx+2,
     1   ny+2,natmi,tmp,'TEXTAN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:natmi-1)
     1   tmp(ix,iy,is)=pfluxa(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:PFLA','s^-1.m^-2',nx+2,
     1   ny+2,natmi,tmp,'TEXTAN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:natmi-1)
     1   tmp(ix,iy,is)=pefluxa(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:PEFA','W.m^-2',nx+2,
     1   ny+2,natmi,tmp,'TEXTAN','')
        forall (ix=-1:nx,iy=-1:ny)
     1   temp(ix,iy)=emiss(ix+1,iy+1,1,1)
        call mds_put_double_rz('\TOP.SNAPSHOT:EMISS','ph.m^-3.s^-1',
     1   nx+2,ny+2,temp,'')
        forall (ix=-1:nx,iy=-1:ny,is=0:natmi-1)
     1   tmp(ix,iy,is)=eneutrad(ix+1,iy+1,is+1,0)
        call mds_put_double_rzns('\TOP.SNAPSHOT:ENEUTRAD','W',nx+2,
     1   ny+2,natmi,tmp,'TEXTAN','')
      endif
#endif
#ifdef B25_EIRENE
      if(nmoli.gt.0) then
        forall (ix=-1:nx,iy=-1:ny,is=0:nmoli-1)
     1   tmp(ix,iy,is)=dmb2(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:DMB2','/m^3',nx+2,
     1   ny+2,nmoli,tmp,'TEXTMN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:nmoli-1)
     1   tmp(ix,iy,is)=tmb2(ix+1,iy+1,is+1,1)/ev
        call mds_put_double_rzns('\TOP.SNAPSHOT:TMB2','eV',nx+2,
     1   ny+2,nmoli,tmp,'TEXTMN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:nmoli-1)
     1   tmp(ix,iy,is)=rfluxm(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:RFLM','s^-1.m^-2',nx+2,
     1   ny+2,nmoli,tmp,'TEXTMN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:nmoli-1)
     1   tmp(ix,iy,is)=refluxm(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:REFM','W.m^-2',nx+2,
     1   ny+2,nmoli,tmp,'TEXTMN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:nmoli-1)
     1   tmp(ix,iy,is)=pfluxm(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:PFLM','s^-1.m^-2',nx+2,
     1   ny+2,nmoli,tmp,'TEXTMN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:nmoli-1)
     1   tmp(ix,iy,is)=pefluxm(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:PEFM','W.m^-2',nx+2,
     1   ny+2,nmoli,tmp,'TEXTMN','')
        forall (ix=-1:nx,iy=-1:ny)
     1   temp(ix,iy)=emissmol(ix+1,iy+1,1,1)
        call mds_put_double_rz('\TOP.SNAPSHOT:EMISSMOL','ph.m^-3.s^-1',
     1   nx+2,ny+2,temp,'')
        forall (ix=-1:nx,iy=-1:ny,is=0:nmoli-1)
     1   tmp(ix,iy,is)=srcml(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:SRCML','/s',nx+2,
     1   ny+2,nmoli,tmp,'TEXTMN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:nmoli-1)
     1   tmp(ix,iy,is)=edissml(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:EDISSML','W',nx+2,
     1   ny+2,nmoli,tmp,'TEXTMN','')
      endif
#endif
#ifdef B25_EIRENE
      if(nioni.gt.0) then
        forall (ix=-1:nx,iy=-1:ny,is=0:nioni-1)
     1   tmp(ix,iy,is)=dib2(ix+1,iy+1,is+1,1)
        call mds_put_double_rzns('\TOP.SNAPSHOT:DIB2','/m^3',nx+2,
     1   ny+2,nioni,tmp,'TEXTIN','')
        forall (ix=-1:nx,iy=-1:ny,is=0:nioni-1)
     1   tmp(ix,iy,is)=tib2(ix+1,iy+1,is+1,1)/ev
        call mds_put_double_rzns('\TOP.SNAPSHOT:TIB2','eV',nx+2,
     1   ny+2,nioni,tmp,'TEXTIN','')
      endif
#endif

      rtmp=(te(jxi,jsep)/qe*vol(jxi,jsep+1)+
     2 te(jxi,jsep+1)/qe*vol(jxi,jsep))/
     3 (vol(jxi,jsep)+vol(jxi,jsep+1))
      call mds_put_double('\top.ident.impsep:te','eV',rtmp)
      dbimpsep_te=rtmp
      rtmp=(ti(jxi,jsep)/qe*vol(jxi,jsep+1)+
     2 ti(jxi,jsep+1)/qe*vol(jxi,jsep))/
     3 (vol(jxi,jsep)+vol(jxi,jsep+1))
      call mds_put_double('\top.ident.impsep:ti','eV',rtmp)
      dbimpsep_ti=rtmp
      rtmp=(ne(jxi,jsep)*vol(jxi,jsep+1)+
     2 ne(jxi,jsep+1)*vol(jxi,jsep))/
     3 (vol(jxi,jsep)+vol(jxi,jsep+1))
      call mds_put_double('\top.ident.impsep:ne','m^-3',rtmp)
      dbimpsep_ne=rtmp
      rtmp=(dna0(jxi,jsep,ismain)*vol(jxi,jsep+1)+
     2 dna0(jxi,jsep+1,ismain)*vol(jxi,jsep))/
     3 (vol(jxi,jsep)+vol(jxi,jsep+1))
      call mds_put_double('\top.ident.impsep:d','m^2/s',rtmp)
      dbimpsep_d=rtmp
      rtmp=(hce0(jxi,jsep)/ne(jxi,jsep)*vol(jxi,jsep+1)+
     2 hce0(jxi,jsep+1)/ne(jxi,jsep+1)*vol(jxi,jsep))/
     3 (vol(jxi,jsep)+vol(jxi,jsep+1))
      call mds_put_double('\top.ident.impsep:kye','m^2/s',rtmp)
      dbimpsep_kye=rtmp
      rtmp=(hcib(jxi,jsep,ismain)/na(jxi,jsep,ismain)*vol(jxi,jsep+1)+
     2 hcib(jxi,jsep+1,ismain)/na(jxi,jsep+1,ismain)*vol(jxi,jsep))/
     3 (vol(jxi,jsep)+vol(jxi,jsep+1))
      call mds_put_double('\top.ident.impsep:kyi','m^2/s',rtmp)
      dbimpsep_kyi=rtmp
      rtmp=(hci0(jxi,jsep)/ni(jxi,jsep,0)*vol(jxi,jsep+1)+
     2 hci0(jxi,jsep+1)/ni(jxi,jsep+1,0)*vol(jxi,jsep))/
     3 (vol(jxi,jsep)+vol(jxi,jsep+1))
      call mds_put_double('\top.ident.impsep:kyi0','m^2/s',rtmp)
      dbimpsep_kyi0=rtmp
      rtmp=cr_y(jxi,jsep+1)
      call mds_put_double('\top.ident.impsep:cr','m',rtmp)
      dbimpsep_cr=rtmp
      rtmp=cz_y(jxi,jsep+1)
      call mds_put_double('\top.ident.impsep:cz','m',rtmp)
      dbimpsep_cz=rtmp

      rtmp=(te(jxa,jsep)/qe*vol(jxa,jsep+1)+
     2 te(jxa,jsep+1)/qe*vol(jxa,jsep))/
     3 (vol(jxa,jsep)+vol(jxa,jsep+1))
      call mds_put_double('\top.ident.ompsep:te','eV',rtmp)
      dbompsep_te=rtmp
      rtmp=(ti(jxa,jsep)/qe*vol(jxa,jsep+1)+
     2 ti(jxa,jsep+1)/qe*vol(jxa,jsep))/
     3 (vol(jxa,jsep)+vol(jxa,jsep+1))
      call mds_put_double('\top.ident.ompsep:ti','eV',rtmp)
      dbompsep_ti=rtmp
      rtmp=(ne(jxa,jsep)*vol(jxa,jsep+1)+
     2 ne(jxa,jsep+1)*vol(jxa,jsep))/
     3 (vol(jxa,jsep)+vol(jxa,jsep+1))
      call mds_put_double('\top.ident.ompsep:ne','eV',rtmp)
      dbompsep_ne=rtmp
      rtmp=(dna0(jxa,jsep,ismain)*vol(jxa,jsep+1)+
     2 dna0(jxa,jsep+1,ismain)*vol(jxa,jsep))/
     3 (vol(jxa,jsep)+vol(jxa,jsep+1))
      call mds_put_double('\top.ident.ompsep:d','m^2/s',rtmp)
      dbompsep_d=rtmp
      rtmp=(hce0(jxa,jsep)/ne(jxa,jsep)*vol(jxa,jsep+1)+
     2 hce0(jxa,jsep+1)/ne(jxa,jsep+1)*vol(jxa,jsep))/
     3 (vol(jxa,jsep)+vol(jxa,jsep+1))
      call mds_put_double('\top.ident.ompsep:kye','m^2/s',rtmp)
      dbompsep_kye=rtmp
      rtmp=(hcib(jxa,jsep,ismain)/na(jxa,jsep,ismain)*vol(jxa,jsep+1)+
     2 hcib(jxa,jsep+1,ismain)/na(jxa,jsep+1,ismain)*vol(jxa,jsep))/
     3 (vol(jxa,jsep)+vol(jxa,jsep+1))
      call mds_put_double('\top.ident.ompsep:kyi','m^2/s',rtmp)
      dbompsep_kyi=rtmp
      rtmp=(hci0(jxa,jsep)/ni(jxa,jsep,0)*vol(jxa,jsep+1)+
     2 hci0(jxa,jsep+1)/ni(jxa,jsep+1,0)*vol(jxa,jsep))/
     3 (vol(jxa,jsep)+vol(jxa,jsep+1))
      call mds_put_double('\top.ident.ompsep:kyi0','m^2/s',rtmp)
      dbompsep_kyi0=rtmp
      rtmp=cr_y(jxa,jsep+1)
      call mds_put_double('\top.ident.ompsep:cr','m',rtmp)
      dbompsep_cr=rtmp
      rtmp=cz_y(jxa,jsep+1)
      call mds_put_double('\top.ident.ompsep:cz','m',rtmp)
      dbompsep_cz=rtmp

      fhesep=0.0_R8
      fhisep=0.0_R8
      fnesep=0.0_R8
      fnisep=0.0_R8
      do iy=-1,ny
        do ix=-1,nx
          if(topix(ix,iy).ge.-1.and.topix(ix,iy).le.nx.and.
     1     topiy(ix,iy).ge.-1.and.topiy(ix,iy).le.ny) then
            if(on_closed_surface(ix,iy).and.
     1         .not.on_closed_surface(topix(ix,iy),topiy(ix,iy))) then
              fhesep=fhesep+fhe(topix(ix,iy),topiy(ix,iy),1)
              fhisep=fhisep+fhi(topix(ix,iy),topiy(ix,iy),1)
              fnesep=fnesep+fne(topix(ix,iy),topiy(ix,iy),1)
              fnisep=fnisep+fni(topix(ix,iy),topiy(ix,iy),1)
            endif
          endif
        enddo
      enddo
      call mds_put_double('\top.ident.sep:fhe','W',fhesep)
      dbsep_fhe=fhesep
      call mds_put_double('\top.ident.sep:fhi','W',fhisep)
      dbsep_fhi=fhisep
      call mds_put_double('\top.ident.sep:fne','#/s',fnesep)
      dbsep_fne=fnesep
      call mds_put_double('\top.ident.sep:fni','#/s',fnisep)
      dbsep_fni=fnisep

      write(*,*) '... b2mds'

#else
      call xerrab('b2mds: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      call dealloc_b2mod_wall
      call dealloc_b2mod_layer
      call dealloc_b2mod_b2plot
      call dealloc_sputter_data
      call dealloc_b2mod_sputter
      call dealloc_b2mod_elements
      call dealloc_b2mod_boundary
      call dealloc_b2mod_neutrals
      call dealloc_b2mod_wall_init
#ifdef B25_EIRENE
      CALL EIRENE_DEALLOC_CTRCEI
      CALL EIRENE_DEALLOC_BRAEIR
      CALL EIRENE_DEALLOC_EIRBRA
      CALL EIRENE_EXTRAB25_CLEANUP
#else
      call dealloc_b2mod_braeir
      call dealloc_b2mod_eirbra
      call dealloc_b2mod_eirdiag
#endif
      call subend()
      return
      end subroutine b2mds
