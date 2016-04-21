      subroutine eirene_mc(nx, ny, ns, ismain, dtim, BoRiS, 
     1 sna0, smo0, she0, shi0, main_call)

c  version : 26.07.2002 23:39

      use b2mod_types
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_geo
      use b2mod_neutrals_namelist
      use b2mod_numerics_namelist
      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_eirene_globals
      use b2mod_tallies
      use b2mod_indirect
      use b2mod_neutr_src_scaling
      use b2mod_feedback
      use b2mod_diag
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_b2cmpb
#ifdef B25_EIRENE
      use b2mod_work
      use eirmod_braeir
      use eirmod_eirbra
      use eirmod_extrab25
#endif
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/BRAEIR.F>
#     include <COUPLE/EIRBRA.F>
#endif
      integer nx,ny,ns,ismain
      real (kind=R8) ::
     * dtim, BoRiS,
     * sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     & she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3)
      logical main_call                                                  !dpc
csw 19oct2012, added save attribute!!!
      logical, save :: my_ltime

      integer, save :: background_new=0                              !srv 14.07.10
      integer ix,iy,is,ib,ih,ifl,istrai,ncall,iatm,imol
      integer pfrregno1, pfrregno2, solregno, isfeedback, ncallfeedback,
     1 nesepm_istra, use_recyceir
      real (kind=R8) ::
     * recyc_hlp,vthe,t0,t1,rf0,uatmp,dua,t2,t3,
     * neutral_sources_rescale, ndes, private_flux_puff, nesepm_sol,
     * recyctmp, nesepm_pfr, volrec_sol, ndes_sol, nepedm_sol
      real (kind=R8), save, allocatable ::
     1 SSNI(:,:), SSMO(:,:), SSEE(:), SSEI(:), SSEICORR(:),
     1 tfl(:,:), volsumnold(:)
      real (kind=R8), save, allocatable ::
     * SMO_TOT(:,:,:), SNI_TOT(:,:,:)

cank-20020724 {
c*** The type of source scaling is determined by the lstrascl parameter
c*** for each stratum, each species:
c***  < 0   - prehistoric treatment (left for compatibility)
c***  = 0   - no source scaling at all
c***  =< nspecies - full treatment, on each sub-iteration
c***  > nspecies  - partial treatment, after Eirene call only
c***
c*** In the Eirene particle balance, we distinguish the following
c*** components for each stratum, each species:
c***
c***  hlp_pi  plasma onto the surfaces
c***  hlp_po  resulting neutrals from the surfaces (recycling)
c***  hlp_ni  neutrals onto the surfaces
c***  hlp_no  reflected neutrals
c***  hlp_g   neutrals from the gas puffing sources
c***  hlp_s   neutrals ionising in the plasma (<0 for recombination)
c***  hlp_t   neutrals absorbed in the "time layer"
c***  tfl     plasma onto the surfaces in B2
c***
c*** For a recycling stratum, first scale all the sources and fluxes by
c*** the ratio of tfl/hlp_pi to have exact match between B2 & Eirene.
c*** Then scale hlp_s, hlp_t, hlp_ni, hlp_no by x ~ 1:
c***
c***    (hlp_s + hlp_t + hlp_ni - hlp_no) * x = hlp_po
c***
c*** (this should only be done for the species scaled with themselves)
c***
c*** For a non-recycling stratum, scale hlp_s & hlp_t by x, and
c*** hlp_ni & hlp_no by 1/x:
c***
c***    (hlp_s + hlp_t) * x + (hlp_ni - hlp_no) / x = hlp_g
c***
c*** selecting the x value positive, closest to 1.
c***
c*** No scaling for the "time layer"
c***
c*** 10.02.2013: If l_neutrad<0, then old (incorrect) treatment of the neutral radiation.
c*** Otherwise: take the radiation explicitly calculated in EIRENE. 

      real (kind=R8) :: s,x1,x2,t00,hlp_st,rec_tol
      real (kind=R8), save, allocatable ::
     1 hlp_pi(:),hlp_po(:),hlp_ni(:),
     2 hlp_no(:),hlp_g(:),hlp_s(:),hlp_t(:),
     1 hlp_scl(:),hlp_scp(:,:),
     2 hlp_stra(:), core_flux(:)
      real (kind=R8), save :: core_surface
      logical hlp_bad_isp, hlp_log, hlp_pr, hlp_tt, hlp_cor
      logical, allocatable, save ::  hlp_pff(:,:)
      integer i,j,k,l,js,ls,kks,n,itmp,icount
      integer, save :: l_ddt
      integer, save, allocatable :: ispc_ref(:), ishigh(:)
      character*13, save, allocatable :: ion_name(:)
      real (kind=R8), allocatable :: sum_wld(:,:)
cank }
      character*24 hlp_frm, hlp_frm1, hlp_frm2, hlp_frm3, hlp_frm4
c      integer lns_hlp
c      parameter (lns_hlp=4)
c      integer lis_hlp(lns_hlp) 
c      data lis_hlp / 110, 111, 112, 49 /
      logical last, saved_volrec_used

      integer output_style, linearisation, ank_mods, save_f31,
     1 repeat_first_call, ionising_core, dpc_fix, eircount,
     2 no_sseicorr, save_f30
      save ncall, output_style, saved_volrec_used,
     * pfrregno1, pfrregno2, solregno, isfeedback,
     * nesepm_sol, nesepm_pfr, ndes, private_flux_puff, nesepm_istra,
     * neutral_sources_rescale, ncallfeedback, ank_mods, save_f31,
     1 repeat_first_call, use_recyceir, ionising_core, dpc_fix,
     2 volrec_sol, ndes_sol, no_sseicorr, nepedm_sol, save_f30
csw 26mar2013
      integer, save :: iscltst=0
csw 26feb2014
      integer, save :: ihist=0, neutr_avg=0, nhist=1, nrelax = 0
      integer, save :: ih_end
      logical, save :: full_history = .false.
      real (kind=R8) :: relax

      real (kind=R8), save, allocatable :: 
     .   sni_av(:,:,:,:),smo_av(:,:,:,:),sei_av(:,:,:),see_av(:,:,:)
      real (kind=R8), allocatable, save ::                        !IYS 28.03.12 {
     &  sni_h(:,:,:,:,:),
     &  smo_h(:,:,:,:,:),
     &  see_h(:,:,:,:), sei_h(:,:,:,:),
     &  pfluxa_h(:,:,:,:,:), rfluxa_h(:,:,:,:,:), 
     &  pfluxm_h(:,:,:,:,:), rfluxm_h(:,:,:,:,:),
! _h is an abbreviation from 'history'
     &  sni_s(:,:,:,:),
     &  smo_s(:,:,:,:),
     &  see_s(:,:,:), sei_s(:,:,:),
     &  pfluxa_s(:,:,:,:), rfluxa_s(:,:,:,:), 
     &  pfluxm_s(:,:,:,:), rfluxm_s(:,:,:,:),
! _s is an abbreviation from 'sum'
     &  dab2_h(:,:,:,:,:), dab2_s(:,:,:,:)

*   ..procedures
      real (kind=R8) :: rratio
      intrinsic min, max, mod, abs, sqrt, log
      external subini, subend, xertst, ipgeti, ipgetr, rratio
      logical lnext
      lnext(is,js) = zamax(is).lt.zamin(js).and.zn(is).eq.zn(js).and.
     &  abs(am(is)-am(js)).lt.1.0e-3_R8
*   ..initialisation
      data ncall/0/, output_style/1/, saved_volrec_used/.false./,
     .  linearisation/1/, solregno/3/,
     1  pfrregno1/0/, pfrregno2/2/, isfeedback/0/, ncallfeedback/0/,
     2  neutral_sources_rescale/1.0_R8/, rec_tol/1.e-4_R8/,
     3  nesepm_sol/0.0_R8/, nesepm_pfr/0.0_R8/, nesepm_istra/-1/,
     4  private_flux_puff/0.0_R8/, ionising_core/0/, dpc_fix/1/,
     5  ndes/0.0_R8/, ank_mods/0/, save_f31/0/, repeat_first_call/1/,
     6  use_recyceir/1/, volrec_sol/0.0_R8/, ndes_sol/0.0_R8/,
     7  no_sseicorr/0/, nepedm_sol/0.0_R8/, save_f30/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini('eirene_mc')
*   ..set internal parameters on first call
      if(ncall.eq.0) then !{
c not sure this is the right place
        allocate (hlp_pi(NATM))
        allocate (hlp_po(NATM))
        allocate (hlp_ni(NATM))
        allocate (hlp_no(NATM))
        allocate (hlp_g(NATM))
        allocate (hlp_s(NATM))
        allocate (hlp_t(NATM))
        allocate (hlp_scl(NATM))
        allocate (hlp_scp(NATM,NSTRAT))
        allocate (hlp_stra(NATM))
        allocate (core_flux(NATM))
        allocate (hlp_pff(NATM,NSTRAT))
        allocate (ion_name(NATM))
        allocate (tfl(NSTRAT,0:ns-1))
        allocate (SSNI(NSTRAT,0:NS-1))
        allocate (SSMO(NSTRAT,0:NS-1)) 
        allocate (SSEE(NSTRAT))
        allocate (SSEI(NSTRAT))  
        allocate (SSEICORR(NSTRAT))
        allocate (SMO_TOT(0:NDXP,0:NDYP,NFL))
        allocate (SNI_TOT(0:NDXP,0:NDYP,NFL))
        allocate (volsumnold(NSTRAT))
        volsumnold=-1.0_R8
        allocate(ispc_ref(NATM))
        allocate(ishigh(NATM))
        allocate(sum_wld(2,2+NATM))

        write(*,*) 'Calling test_couple'
        call test_couple(ndx,ndy,nstra,nfl,natm,nmol,nion,nlimps)
        write(*,*) 'Called test_couple'
        call ipgeti ('eirene_mc_output_style', output_style)
        call ipgeti ('eirene_repeat_first_call', repeat_first_call)
        call ipgeti ('eirene_mc_linearisation', linearisation)
        if (linearisation.eq.1)   ! for backward compatibility
     >   call ipgeti ('eirene_mc_linearization', linearisation)
        call xertst (linearisation.eq.1 .or. linearisation.eq.0,
     .   'faulty input parameter linearisation')
        call ipgetr ('b2stbc_ndes', ndes)
        call xertst (0.0_R8.le.ndes, 'faulty internal parameter ndes')
        call ipgetr ('b2stbc_nesepm_sol', nesepm_sol)
        call xertst (0.0_R8.le.nesepm_sol,
     1   'faulty internal parameter nesepm_sol')
        call ipgetr ('b2stbc_nepedm_sol', nepedm_sol)
        call xertst (0.0_R8.le.nepedm_sol,
     1   'faulty internal parameter nepedm_sol')
        call ipgetr ('b2stbc_nesepm_pfr', nesepm_pfr)
        call xertst (0.0_R8.le.nesepm_pfr,
     1   'faulty internal parameter nesepm_pfr')
        call ipgetr ('b2stbc_volrec_sol', volrec_sol)
        call xertst (0.0_R8.le.volrec_sol,
     1   'faulty internal parameter volrec_sol')
        call ipgetr ('b2stbc_ndes_sol', ndes_sol)
        call xertst (0.0_R8.le.ndes_sol,
     1   'faulty internal parameter ndes_sol')
        call ipgeti ('eirene_nesepm_istra', nesepm_istra)
        if(nesepm_istra.gt.0) then
          call xertst (nesepm_istra.le.nstrat,
     1     'faulty internal parameter nesepm_istra')
          call xertst (crcstra(nesepm_istra).eq.'C',
     1     'Stratum nesepm_istra is not declared as a puff stratum!')
        else
          do istrai = 1, nstrat
            if(nesepm_istra.le.0.and.CRCSTRA(istrai).EQ.'C') then
              nesepm_istra=istrai
            endif
          enddo
        endif
        if (ndes.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8.or.
     &    private_flux_puff.gt.0.0_R8.or.nesepm_sol.gt.0.0_R8.or.
     *    volrec_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8.or.
     &    nepedm_sol.gt.0.0_R8) then
          call xertst (nesepm_istra.gt.0,
     &     'No "C" stratum for feedback scheme specified !')
          icount = 0
          do istrai = 1, nstrat
            if (crcstra(istrai).eq.'C' .and.
     .          species_start(istrai).le.isfeedback .and.
     .          species_end(istrai).ge.isfeedback) icount = icount+1
          end do
          call xertst (icount.eq.1 .or. nesepm_istra.gt.0,
     .     'Multiple definitions of puff strata, use nesepm_istra') 
        endif
        call ipgetr ('b2stbc_private_flux_puff', private_flux_puff)
        call xertst (0.0_R8.le.private_flux_puff,
     1   'faulty internal parameter private_flux_puff')
        call ipgeti ('b2stbc_pfrregno1', pfrregno1)
        call ipgeti ('b2stbc_pfrregno2', pfrregno2)
        call ipgeti ('b2stbc_solregno', solregno)
        call ipgeti ('b2stbc_isfeedback', isfeedback)
        call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1   'faulty internal parameter isfeedback')
        if (ndes.gt.0.0_R8.or.private_flux_puff.gt.0.0_R8
     1     .or.nesepm_pfr.gt.0.0_R8) then !{
          call xertst (cbirso.le.pfrregno1.and.
     1     pfrregno1.le.cbirso+cbnrso-1,
     2     'faulty internal parameter pfrregno1')
          call xertst (cbirso.le.pfrregno2.and.
     1     pfrregno2.le.cbirso+cbnrso-1,
     2     'faulty internal parameter pfrregno2')
        endif !}
        if (nesepm_sol.gt.0.0_R8.or.volrec_sol.gt.0.0_R8.or.
     1      nepedm_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8) then
          call xertst (cbirno.le.solregno.and.
     1     solregno.le.cbirno+cbnrno-1,
     2     'faulty internal parameter solregno')
        endif
        call ipgetr ('b2mndr_rescale_neutrals_sources',
     *       neutral_sources_rescale)
        call xertst (0.0_R8.le.neutral_sources_rescale,
     1   'faulty input neutral_sources_rescale')
        call ipgeti ('b2stbc_ncallfeedback', ncallfeedback)
        call ipgeti ('eirene_ionising_core', ionising_core)
        if (ionising_core.eq.0)     ! for backward compatibility
     >   call ipgeti ('eirene_ionizing_core', ionising_core)
        call xertst (ionising_core.gt.0.or.abs(ionising_core).le.nbc,
     1   'faulty input ionising_core, abs(ionising_core) > NBC')
        call ipgeti ('eirene_ank_mods', ank_mods)
        call ipgeti ('eirene_dpc_fix', dpc_fix)
        call ipgeti ('eirene_savef30', save_f30)
        call ipgeti ('eirene_savef31', save_f31)
        call ipgeti ('eirene_use_recyceir', use_recyceir)
csw 28jan2011, call to eirene_mc_init not used anymore, use new modules
c        call eirene_mc_init
        call eirene_extrab25_alloc_mods
!pb 13dec2013 eirpbls removed - not needed
!pb        call eirene_extrab25_eirpbls_init(1,
!pb     .                       natm,nmol,nion,npls,nstra,0)
!pb        call eirene_extrab25_eirpbls_init(2,0,0,0,0,nstra,npls)
csw 22jul2011
        call eirene_extrab25_iniusr_init(n_spcsrf,l_spcsrf,i_spcsrf,
     .           j_spcsrf,sps_sgrp,sps_absr,sps_trno,sps_trni,
     .           sps_mtri,sps_tmpr,sps_spph,sps_spch,
     .           sps_mtrl,sps_id)
#ifndef NO_B2_CHEM_SPUT
        call eirene_extrab25_braeir_init(chemical_sputter_yield)
#endif
csw 26mar2013
        call ipgeti ('eirene_iscltst',iscltst)

csw 26feb2013
        call ipgeti ('eirene_neutr_avg',neutr_avg)
        call ipgeti ('b2stbr_eir_src_nhist',nhist)
        call ipgeti ('eirene_underrelax',nrelax)
        call xertst (.not.(
     .    (neutr_avg.gt.1 .and. (nhist.gt.1 .or. nrelax.gt.1)) .or.
     .    (nhist.gt.1 .and. (neutr_avg.gt.1 .or. nrelax.gt.1)) .or.
     .    (nrelax.gt.1 .and. (neutr_avg.gt.1 .or. nhist.gt.1))),
     .   'Cannot have simultaneous multiple Eirene averaging schemes!')
        if(nhist.gt.1) call alloc_hist_arrays(nhist)
        if(nrelax.gt.1 .or. nhist.gt.1) call alloc_sum_arrays
        if(lhalpha) then
          itmp=1
        else
          itmp=0
        endif
        call ipgeti('eirene_lhalpha',itmp)
        lhalpha=itmp.ne.0
        write(*,*) 'lhalpha = ',lhalpha
        if(lvib) then
          itmp=1
        else
          itmp=0
        endif
        call ipgeti('eirene_lvib',itmp)
        lvib=itmp.ne.0
        write(*,*) 'lvib = ',lvib
        call xertst(.not.(ank_mods.ne.0 .and. ionising_core.ne.0),
     1   'contradictory inputs ank_mods and ionising_core!')
        if(eirene_mod.le.0) eirene_mod=1
        call ipgeti ('eirene_no_sseicorr', no_sseicorr)
        call ipgeti ('eirene_background_new', background_new)        !srv 14.07.10
      end if !}

      if(ncall.eq.0) then
        call write_f30(nx, ny)
      endif
      if(ncall.eq.0.or.(main_call.and.mod(ncall,eirene_mod).eq.0)) then
        if(background_new.eq.0) call copy_background(nx,ny,ns,BoRiS)  !srv 14.07.10
        if(background_new.ne.0) call copy_background_(nx,ny,ns,BoRiS) !srv 14.07.10
        if(ncall.eq.0.or.save_f31.ne.0) then
          call write_f31(nx, ny)
        endif
      endif

      if(ncall.eq.0) then !{
        call alloc_b2mod_neutrals_save(nstraid)
        call read_b2mod_neutrals_save(99,'    ',saved_volrec_used)
        call xertst(nstrat.le.nstraid,
     &   'Increase DEF_NSTRA in DIMENSIONS.F')
        if(maxval(saved_volrec).gt.0.0_R8) then !{
          volsumnold(1:nstrat) = saved_volrec(1:nstrat)
        else !}{
          do istrai = 1, nstrat !{
            if (CRCSTRA(istrai).EQ.'V') then !{
              if(volsumnold(istrai).lt.0.0_R8) then
                volsumnold(istrai)=volrecstart(istrai)
                saved_volrec_used =.true.
              endif
            else !}{
              volsumnold(istrai)=0.0_R8
            end if !}
          end do !}
        end if !}

        if(ank_mods.ne.0 .or. ionising_core.ne.0) then
c*** Check the b2espcr array
          hlp_bad_isp=.false.
          ispc_ref=-1
          do i=0,ns-1 !{
cwd exclude fluid neutrals from this check: their recycling coefficient does not
c   necessarily have to match that of the ions, since they are not treated in the
c   same Eirene stratum
            if (is_neutral(i)) cycle
            j=b2espcr(i)
            if(ispc_ref(j).lt.0) then !{
              ispc_ref(j)=i
              do k=1,nstrat !{
                if(use_recyceir.gt.0) then
                  recyctmp=recyceir(k)
                else
                  recyctmp=recyc(i,k)
                endif
                hlp_scp(j,k)=recyctmp
              end do !}
            else !}{
              hlp_tt=.false.
              do k=1,nstrat !{
                if(use_recyceir.gt.0) then
                  recyctmp=recyceir(k)
                else
                  recyctmp=recyc(i,k)
                endif
                hlp_tt=hlp_tt.or.
     1           abs(hlp_scp(j,k)-recyctmp).gt.rec_tol
              end do !}
              if(hlp_tt) then !{
                write (0,*)
     1           'eirene_mc: different values of recyceir ',
     ,                                  'species,is= ',j,i
                hlp_bad_isp=.true.
              end if !}
            end if !}
          end do !}
          call xertst(.not.hlp_bad_isp,
     1     'eirene_mc: error in recyceir specification')

csw 08feb2011 l_ddt definition moved down after eirene call when it's
c             defined
c c*** "time layer"
c          l_ddt=0
c          if(crcstra(nstrat).eq.'T') l_ddt=nnlimi+nnstsi

        end if !}

        hlp_tt=.false.
      endif

      if (ncall.eq.0 .and. ionising_core.ne.0) then
cxpb compute core boundary surface area
        core_surface = 0.0_R8
        do iy = -1, ny
          do ix = -1, nx
            if (bottomiy(ix,iy).eq.-2 .and. mod(region(ix,iy,0),4).eq.1)
     >       core_surface = core_surface + 
     .        gs(topix(ix,iy),topiy(ix,iy),1)
          enddo
        enddo
cxpb compute ishigh and ion_name arrays
        ishigh = -1
        ion_name = ' '
        do is = 0, ns-1
          j = b2espcr(is)
          if (ishigh(j).eq.-1) then
            do js = is+1, ns-1
              if (lnext(is,js)) ishigh(j) = js
            enddo
            if (ishigh(j).eq.-1) ishigh(j) = is
            call species (ishigh(j), ion_name(j), .false.)
          endif
        enddo
      endif

*  -- COMPUTE TOTAL FLUXES TO RECYCLING BOUNDARIES --

c*** tflux(i)  : source strength for stratum i as given to Eirene
c*** tfl(i,j)  : the same for species j only
c*** flux_scale: conversion factor Eirene-B2
c*** recyc     : fudge factors for the sources in B2

      tfl=0.0_R8
      tflux=0.0_R8
      flux_scale=0.0_R8
      do istrai = 1, nstrat !{
        flux_scale(istrai)=1.0_R8
        IF(crcstra(istrai).eq.'W') then !{
          do is = species_start(istrai),species_end(istrai) !{
            if(.not.is_neutral(is)) then !{
              do iy=1,rc_list_size(istrai) !{
                tflux(istrai)=tflux(istrai)-
     1      fna(rightix(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),
     2          rightiy(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),0,is)
                tfl(istrai,b2espcr(is))=tfl(istrai,b2espcr(is))-
     1      fna(rightix(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),
     2          rightiy(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),0,is)
              end do !}
            end if !}
          end do !}
        else if(crcstra(istrai).eq.'E') then !}{
          do is = species_start(istrai),species_end(istrai) !{
            if(.not.is_neutral(is)) then !{
              do iy=1,rc_list_size(istrai) !{
                tflux(istrai)=tflux(istrai)+
     .           fna(rc_list_x(iy,istrai),rc_list_y(iy,istrai),0,is)
                tfl(istrai,b2espcr(is))=tfl(istrai,b2espcr(is))+
     1           fna(rc_list_x(iy,istrai),rc_list_y(iy,istrai),0,is)
              end do !}
            end if !}
          end do !}
        else if(crcstra(istrai).eq.'S') then !}{
          do is = species_start(istrai),species_end(istrai) !{
            if(.not.is_neutral(is)) then !{
              do ix=1,rc_list_size(istrai) !{
                tflux(istrai)=tflux(istrai)-
     1      fna(topix(rc_list_x(ix,istrai),rc_list_y(ix,istrai)),
     2          topiy(rc_list_x(ix,istrai),rc_list_y(ix,istrai)),1,is)
                tfl(istrai,b2espcr(is))=tfl(istrai,b2espcr(is))-
     1      fna(topix(rc_list_x(ix,istrai),rc_list_y(ix,istrai)),
     2          topiy(rc_list_x(ix,istrai),rc_list_y(ix,istrai)),1,is)
              end do !}
            end if !}
          end do !}
        else if(crcstra(istrai).eq.'N') then !}{
          do is = species_start(istrai),species_end(istrai) !{
            if(.not.is_neutral(is)) then !{
              do ix=1,rc_list_size(istrai) !{
                tflux(istrai)=tflux(istrai)+
     .           fna(rc_list_x(ix,istrai),rc_list_y(ix,istrai),1,is)
                tfl(istrai,b2espcr(is))=tfl(istrai,b2espcr(is))+
     1           fna(rc_list_x(ix,istrai),rc_list_y(ix,istrai),1,is)
              end do !}
            end if !}
          end do !}
        else if(crcstra(istrai).eq.'C') then !}{
          if (ncall.ge.ncallfeedback.and.istrai.eq.nesepm_istra.and.
     &     (ndes.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8.or.
     &      private_flux_puff.gt.0.0_R8)) then !{
            if (pfrregno1.ne.pfrregno2) then !{
              tflux(istrai)=(cbsna(0,isfeedback,pfrregno1)+
     1         cbsna(0,isfeedback,pfrregno2))/neutral_sources_rescale
            else !}{
              tflux(istrai)=
     1         cbsna(0,isfeedback,pfrregno1)/neutral_sources_rescale
            end if !}
            if(ncall.eq.0) then
              if(nesepm_pfr.gt.0.0_R8) then
                write(*,*) '"nesepm_pfr" feedback on stratum ',istrai
              else if(ndes.gt.0.0_R8) then
                write(*,*) '"ndes" feedback on stratum ',istrai
              else if(private_flux_puff.gt.0.0_R8) then
                write(*,*) '"private_flux_puff" feedback on stratum ',
     1           istrai
              endif
            endif
          else if (ncall.ge.ncallfeedback.and.nesepm_sol.gt.0.0_R8.and.
     1     istrai.eq.nesepm_istra) then !}{
            TFLUX(istrai)=
     1       (cbsna(0,isfeedback,solregno))/neutral_sources_rescale
            if(ncall.eq.0)
     1       write(*,*) '"nesepm_sol" feedback on stratum ',istrai
          else if (ncall.ge.ncallfeedback.and.volrec_sol.gt.0.0_R8.and.
     1     istrai.eq.nesepm_istra) then !}{
            TFLUX(istrai)=
     1       (cbsna(0,isfeedback,solregno))/neutral_sources_rescale
            if(ncall.eq.0)
     1       write(*,*) '"volrec_sol" feedback on stratum ',istrai
          else if (ncall.ge.ncallfeedback.and.ndes_sol.gt.0.0_R8.and.
     1     istrai.eq.nesepm_istra) then !}{
            TFLUX(istrai)=
     1       (cbsna(0,isfeedback,solregno))/neutral_sources_rescale
            if(ncall.eq.0)
     1       write(*,*) '"ndes_sol" feedback on stratum ',istrai
          else if (ncall.ge.ncallfeedback.and.nepedm_sol.gt.0.0_R8.and.
     1     istrai.eq.nesepm_istra) then !}{
            TFLUX(istrai)=
     1       (cbsna(0,isfeedback,solregno))/neutral_sources_rescale
            if(ncall.eq.0)
     1       write(*,*) '"nepedm_sol" feedback on stratum ',istrai
          else !}{
            tflux(istrai)=userfluxparm(istrai,1)
          end if !}
          do is = species_start(istrai),species_end(istrai) !{
            tfl(istrai,is)=1.0_R8
          end do !}
          if(tflux(istrai).gt.0.0_R8) then !{
            flux_scale(istrai)=1.0_R8/(tflux(istrai)*qe)
          else !}{
            flux_scale(istrai)=0.0_R8
          end if !}
        else if(crcstra(istrai).eq.'V') then !}{
          tflux(istrai)=1.0_R8/qe
          do is = species_start(istrai),species_end(istrai) !{
            tfl(istrai,is)=1.0_R8
          end do !}
        else if(crcstra(istrai).eq.'T') then !}{
          tflux(istrai)=1.0_R8/qe
          do is = species_start(istrai),species_end(istrai) !{
            tfl(istrai,is)=1.0_R8
          end do !}
c        else if(crcstra(istrai).eq.'U') then !}{
cc         tflux(istrai)=userfluxparm(istrai,1)
c          write(*,*) 'Option "U" not yet implemented'
c          tflux(istrai)=0.0_R8
c          do is = species_start(istrai),species_end(istrai) !{
c            tfl(istrai,is)=1.0_R8
c          end do !}
c          if(tflux(istrai).gt.0.0_R8) then !{
cc warning: the following looks strange but is apparently correct
c             flux_scale(istrai)=1.0_R8/(tflux(istrai)*qe)
c          else !}{
c             flux_scale(istrai)=0.0_R8
c          end if !}
        else !}{
          WRITE(*,*) 'OPTION ',CRCSTRA(istrai),' NOT AVAILABLE'
        end if !}
      end do !}

      write(hlp_frm,'(a,i2,a)') '(a,1p,',nstrat,'(1x,1e14.7))'
      if(main_call.or.ncall.eq.0.or.output_style.gt.1) then !{
        write(*,hlp_frm) 'tflux: ',tflux
        write(*,hlp_frm) 'flux_scale: ',flux_scale
      end if !}

      last=.false.

cxpb !!!!!!!!!!!!!!!!!!!!!!!!!!! EIRENE CALL HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      eircount=0
      if(main_call .and. mod(ncall,eirene_mod).eq.0) then !{

csw 25jan2011 LTIME=.false. in eirsrt call
        my_ltime=.false.

        call eirene_eirsrt(last,my_ltime,eirene_step_dt,tflux,
     *   b2brem,b2rad,b2qie(0),b2vdp,
     1   eirene_step_cpu)
        eircount=1
        if(ncall.eq.0) then 
          write(*,*) 'Calling post_couple'
          call post_couple(nx,ny,ns,nstrai)
          write(*,*) 'Called post_couple'
        endif
        if(ncall.eq.0.and.repeat_first_call.gt.0) then
          do while (eircount.le.repeat_first_call)
           if(background_new.eq.0) call copy_background(nx,ny,ns,BoRiS)  !srv 14.07.10
           if(background_new.ne.0) call copy_background_(nx,ny,ns,BoRiS) !srv 14.07.10

csw 25jan2011 LTIME=.false. in eirsrt call
           my_ltime=.false.
           call eirene_eirsrt(last,my_ltime,eirene_step_dt,tflux,
     *       b2brem,b2rad,b2qie(0),b2vdp,
     1       eirene_step_cpu)
           eircount=eircount+1
          enddo
        endif
      elseif(.not.main_call) then !}{
        if (ncall.lt.3)
     1   write(*,*) 'Skipping call to EIRENE because not main_call'
      elseif(mod(ncall,eirene_mod).ne.0) then !}{
        write(*,*) 'Skipping call to EIRENE because of mod condition'
      end if !}
cxpb !!!!!!!!!!!!!!!!!!!!!!!!!!! EIRENE CALL ENDED !!!!!!!!!!!!!!!!!!!!!!!!!!!!

csw 08feb2011 put l_ddt definition here since nnlimi(nlim) is known only after eirene call
      if(ank_mods.ne.0 .or. ionising_core.eq.1) then
c*** "time layer"
        l_ddt=0
        if(crcstra(nstrai).eq.'T') then
          l_ddt=nnlimi+nnstsi
csw 2011 in case of no time layer in eirene
          if(.not.my_ltime) then
csw 2011
c           l_ddt=l_ddt+1
            l_ddt=0
          endif
        endif
      end if !}

      if(main_call.or.ncall.eq.0.or.output_style.gt.1) then !{
        write(*,hlp_frm) 'tflux: ',tflux
        write(*,hlp_frm) 'flux_scale: ',flux_scale
      end if !}

cxpb  Additional diagnostic for comparison with SOLPS4.0
      if (ncall.eq.0) call b2trci(nx,ny,ns,ismain,dtim,BoRiS)
      if (main_call.and.mod(ncall,eirene_mod).eq.0
     & .and.ank_tracing.ne.0.and.output_style.ge.2) then
        write(*,*) 'Eirene output in eirene_mc :'
        sum_wld(1:2,1:2+nspecies)=0.0_R8
        do istrai = 1, nstrat
          sum_wld(1,1)=sum_wld(1,1)+wldnek(nnlimi+2,istrai)
          sum_wld(2,1)=sum_wld(2,1)+wldnek(nnlimi+3,istrai)
          sum_wld(1,2)=sum_wld(1,2)+wldnep(nnlimi+2,istrai)
          sum_wld(2,2)=sum_wld(2,2)+wldnep(nnlimi+3,istrai)
          do is = 1, nspecies
            sum_wld(1,2+is)=sum_wld(1,2+is)+wldna(nnlimi+2,is,istrai)
            sum_wld(2,2+is)=sum_wld(2,2+is)+wldna(nnlimi+3,is,istrai)
          end do
        end do
        write(*,*) 'Summed inner target data : '
        write(*,*) 
     .   'wldnek, wldnep, wldna = ', sum_wld(1,1:2+nspecies)
        write(*,*) 'Summed outer target data : '
        write(*,*) 
     .   'wldnek, wldnep, wldna = ', sum_wld(2,1:2+nspecies)
      endif
cxpb

cxpb  Compute recombination sources rescaling factor if necessary
      do istrai = 1, nstrat !{
        if(volsumnold(istrai).gt.0.0_R8) then !{
          if(output_style.gt.1.or.main_call) 
     .     write(*,'(a,i2,1p,2(1x,1e14.7))')
     .      'stratum, volsumnold, abs(volsumn) ',
     .      istrai, volsumnold(istrai), abs(volsumn(istrai))/qe
        end if !}
        if(volsumnold(istrai).gt.0.0_R8.and.volrecinc.gt.0.0_R8) then !{
          if(abs(volsumn(istrai)/qe).gt.
     &           volsumnold(istrai)*volrecinc) then !{
            if(main_call) then !{
              flux_scale(istrai)=volsumnold(istrai)*volrecinc/
     1         abs(volsumn(istrai)/qe)
              if(output_style.ge.1) write(*,'(a,i2,a,1p,g12.5)')
     .         'Rescaling volume recombination stratum ',
     1         istrai,' by ',flux_scale(istrai)
            else !}{
              flux_scale(istrai)=volsumnold(istrai)/
     .         abs(volsumn(istrai)/qe)
            end if !}
          end if !}
          volsumnold(istrai)=(1.0_R8-volrecwt)*volsumnold(istrai)+
     1     volrecwt*abs(volsumn(istrai)/qe)*flux_scale(istrai)
          saved_volrec(istrai) = volsumnold(istrai)
        end if !}
      end do !}
cxpb

      if(main_call.and.mod(ncall,eirene_mod).eq.0.and.ank_mods.ne.0)
     > then
        if(ncall.eq.0) then !{
          hlp_tt=.false.
          hlp_bad_isp=.false.
          do i=0,ns-1 !{
            hlp_tt=b2espcr(i).lt.0 .or. b2espcr(i).gt.nnatmi
            hlp_bad_isp= hlp_bad_isp .or. hlp_tt
            if(hlp_tt) write (0,*) 'is,b2espcr,nnatmi: ',
     ,         i,b2espcr(i),nnatmi
          end do !}
          call xertst (.not.hlp_bad_isp,
     1       'eirene_mc: wrong value found in b2espcr')
        end if !}

        if(dbg_eir_mc.gt.9 .and. 
     &       main_call .and. mod(ncall,eirene_mod).eq.0) then !{
          print *
          print *,'%%% Eirene data in B2 %%%',nlimps
          print '(a,5(i3,1x))','nnatmi,nnmoli,nnplsi,nnstsi,nnlimi = ',
     ,       nnatmi,nnmoli,nnplsi,nnstsi,nnlimi
          write(hlp_frm1,'(a,i2,a)')
     ,     '(1x,',2*nnatmi+2*nnmoli,'(6x,a2,i3.2))'
          write(hlp_frm2,'(a,i2,a)')
     ,     '(1p,i4,',2*nnatmi+2*nnmoli,'e11.4)'
          write(hlp_frm3,'(a,i2,a)')
     ,     '(1p,i4,',nnplsi,'e11.4)'
          write(hlp_frm4,'(a,i2,a)')
     ,     '(1p,i4,',1+nnatmi+nnmoli,'e11.4)'
          do k=0,nstrat !{
            print *
            print *,'wldna, wldnm, wldra, wldrm:',k
            write (*,hlp_frm1)
     ,         ('na',i,i=1,nnatmi),('nm',i,i=1,nnmoli),
     ,         ('ra',i,i=1,nnatmi),('rm',i,i=1,nnmoli)
            do j=1,nnlimi+nnstsi !{
              if(j.le.nnlimi .or. j.gt.nnlimi) then !{
                hlp_pr=.false.
                do i=1,nnatmi !{
                  hlp_pr= hlp_pr .or. wldna(j,i,k).ne.0.0_R8
     .               .or. wldra(j,i,k).ne.0.0_R8
                end do !}
                do i=1,nnmoli !{
                  hlp_pr= hlp_pr .or. wldnm(j,i,k).ne.0.0_R8
     .               .or. wldrm(j,i,k).ne.0.0_R8
                end do !}
                if(hlp_pr) then !{
                  write(*,hlp_frm2) j,
     ,               (wldna(j,i,k),i=1,nnatmi),
     ,               (wldnm(j,i,k),i=1,nnmoli),
     ,               (wldra(j,i,k),i=1,nnatmi),
     ,               (wldrm(j,i,k),i=1,nnmoli)
                end if !}
              end if !}
            end do !}

            print *
            print *,'wldpp:',k
            do j=nnlimi+1,nnlimi+nnstsi !{
              hlp_pr=.false.
              do i=1,nnplsi !{
                hlp_pr = hlp_pr .or. wldpp(j,i,k).ne.0.0_R8
              end do !}
              if (hlp_pr)
     >         write(*,hlp_frm3) j,(wldpp(j,i,k),i=1,nnplsi)
            end do !}
            print *
            print *,'wldpeb, wldpa, wldpm:',k
            do j=nnlimi+1,nnlimi+nnstsi !{
              hlp_pr=wldpeb(j,k).ne.0.0_R8
              do i=1,nnatmi !{
                hlp_pr= hlp_pr .or. wldpa(j,i,k).ne.0.0_R8
              end do !}
              do i=1,nnmoli !{
                hlp_pr= hlp_pr .or. wldpm(j,i,k).ne.0.0_R8
              end do !}
              if (hlp_pr) write(*,hlp_frm4) j,wldpeb(j,k),
     ,           (wldpa(j,i,k),i=1,nnatmi),
     ,           (wldpm(j,i,k),i=1,nnmoli)
            end do !}
          end do !}
          print *,'%%%'
        end if !}
        if(dbg_eir_mc.gt.0 .and.
     &       main_call .and. mod(ncall,eirene_mod).eq.0) then !{
c*** Print the particle fluxes for each stratum and species
          print '(a)',
     ,     '%%%=============== total fluxes by stratum and '//
     ,     'species before correction ===============%%%'
          print '(a25,20(4x,a5,4x))',' Stratum Species lstrascl',
     ,       ' tfl ','wldpi','wldpo','wldni','wldno',
     ,       'Srcpl','puff ','d/dt ','crsrc','error'
          do istrai=1,nstrat !{
            do ls=1,nnatmi !{
              hlp_pi(ls)=0.0_R8
              hlp_po(ls)=0.0_R8
              hlp_ni(ls)=0.0_R8
              hlp_no(ls)=0.0_R8
              hlp_g(ls) =0.0_R8
              hlp_s(ls) =0.0_R8
              hlp_t(ls) =0.0_R8
            end do !}
            do j=1,nnlimi+nnstsi !{
              if(j.gt.nnlimi .and. j.ne.l_ddt) then !{
                do is=0,ns-1 !{
                  if(.not.is_neutral(is)) hlp_pi(b2espcr(is))=
     +             hlp_pi(b2espcr(is))+wldpp(j,fluids_list(is),istrai)
                end do !}
                do i=1,nnatmi !{
                  hlp_po(i)=hlp_po(i)+wldpa(j,i,istrai)
                  do is=1,nnmoli !{
                    s=mlcmp(i,is)
                    if(s.gt.0.0_R8)
     -                 hlp_po(i)=hlp_po(i)+wldpm(j,is,istrai)*s
                  end do !}
                end do !}
              end if !}
              hlp_log=j.le.nnlimi
              if(hlp_log .or. (j.gt.nnlimi .and. j.ne.l_ddt)) then !{
                if(.not.hlp_log) then !{
                  do i=1,ncns !{
                    hlp_log=hlp_log .or. j.eq.lcns(i)+nnlimi
                  end do !}
                  do i=1,ntns !{
                    hlp_log=hlp_log .or. j.eq.ltns(i)+nnlimi
                  end do !}
                end if !}

                if(hlp_log) then !{
                  do i=1,nnatmi !{
                    hlp_ni(i)=hlp_ni(i)+wldna(j,i,istrai)
                    hlp_no(i)=hlp_no(i)+wldra(j,i,istrai)
                    do is=1,nnmoli !{
                      s=mlcmp(i,is)
                      hlp_ni(i)=hlp_ni(i)+wldnm(j,is,istrai)*s
                      hlp_no(i)=hlp_no(i)+wldrm(j,is,istrai)*s
                    end do !}
                  end do !}
                end if !}
              else if(j.eq.l_ddt) then !}{
                do i=1,nnatmi !{
                  hlp_t(i)=hlp_t(i)+wldna(j,i,istrai)
                  do is=1,nnmoli !{
                    s=mlcmp(i,is)
                    hlp_t(i)=hlp_t(i)+wldnm(j,is,istrai)*s
                  end do !}
                end do !}
              end if !}
            end do !}

            do ix=1,nx !{
              do iy=1,ny !{
                if(resignore(ix-1,iy-1,1).ne.0) then !{
                  do is=0,ns-1 !{
                    if(.not.is_neutral(is)) then !{
                      hlp_s(b2espcr(is))=hlp_s(b2espcr(is))+
     +                   sni(ix,iy,fluids_list(is),istrai)
                    end if !}
                  end do !}
                end if !}
              end do !}
            end do !}

            do ls=1,nnatmi !{
              kks=ispc_ref(ls)
              hlp_s(ls)=hlp_s(ls)*tflux(istrai)*flux_scale(istrai)
              if(hlp_pi(ls).gt.0.0_R8) then !{
                t0=tfl(istrai,ls)/hlp_pi(ls)
              else !}{
                t0=1.0_R8
              end if !}
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(istrai)
              else
                recyctmp=recyc(kks,istrai)
              endif
              if(crcstra(istrai).eq.'C')
     -           hlp_g(ls)=tflux(istrai)*
     *           gpfc(ls,istrai)*recyctmp
              s=hlp_po(ls)+hlp_g(ls)+hlp_no(ls)-
     -           (hlp_ni(ls)+hlp_s(ls)+hlp_t(ls))
              print '(i6,2i8,3x,1p,10e13.5)',istrai,ls,
     ,           lstrascl(istrai,ls),tfl(istrai,ls),
     ,           hlp_pi(ls),hlp_po(ls),hlp_ni(ls),hlp_no(ls),
     ,           hlp_s(ls),hlp_g(ls),hlp_t(ls),t0,s
            end do !}
          end do !}
          print '(a)',
     ,     '%%%============================================'//
     ,     '=========================================%%%'
        end if !}
      endif

c*** Re-normalisation of the sources and fluxes in order to enforce
c*** the particle balance for recycling and recombining species

csw 23nov2012
      hlp_tt=.false.

csw 2012 fixed rescaling procedure
c        (otherwise rescaling procedure does not correctly work with internal iterations nstg2>1)
c        this is phase 1 of rescaling process
c
      if(main_call.and.mod(ncall,eirene_mod).eq.0
     .            .and.ank_mods.ne.0) then
        write(hlp_frm1,'(a,i2,a)') '(a,1p,i3,',nnatmi,'e13.5)'
        do istrai=1,nstrat !{
          if((crcstra(istrai).eq.'W' .or.
     .        crcstra(istrai).eq.'E' .or.
     .        crcstra(istrai).eq.'N' .or.
     .        crcstra(istrai).eq.'S')) then

c*** Recycling stratum
c*** First, correct the balance between B2 and Eirene sources for the
c*** species to be scaled with themselves (recycling species)

            do ls=1,nnatmi !{
              hlp_pi(ls)=0.0_R8
              hlp_po(ls)=0.0_R8
              hlp_ni(ls)=0.0_R8
              hlp_no(ls)=0.0_R8
              hlp_g(ls) =0.0_R8
              hlp_s(ls) =0.0_R8
              hlp_t(ls) =0.0_R8
            end do !}
            n=nnstsi
            if(l_ddt.gt.0) n=n-1
            do i=1,n !{
              j=nnlimi+i
              do is=0,ns-1 !{
                if(.not.is_neutral(is)) then !{
                  ls=b2espcr(is)
                  if(lstrascl(istrai,ls).eq.ls) then !{
                    hlp_pi(ls)=hlp_pi(ls)+
     .               wldpp(j,fluids_list(is),istrai)
                  end if !}
                end if !}
              end do !}
            end do !}
            do ls=1,nnatmi !{
              if(hlp_pi(ls).gt.0.0_R8) then !{
                hlp_scl(ls)=tfl(istrai,ls)/hlp_pi(ls)
              else !}{
                hlp_scl(ls)=1.0_R8
              end if !}
            end do !}
c!###{
            if(dbg_eir_mc.gt.0) then
              write(*,hlp_frm1) ' @@@ istrai,pi: ',
     ,                            istrai,(hlp_pi(ls),ls=1,nnatmi)
              write(*,hlp_frm1) ' @@@       tfl: ',
     ,                        istrai,(tfl(istrai,ls),ls=1,nnatmi)
              write(*,hlp_frm1) ' @@@       scl: ',
     ,                           istrai,(hlp_scl(ls),ls=1,nnatmi)
            endif
c!###}

            do ix=1,nx !{
              do iy=1,ny !{
                do is=0,ns-1 !{
                  if(.not.is_neutral(is)) then !{
                    ls=b2espcr(is)
                    if(lstrascl(istrai,ls).gt.0) then !{
                      sni(ix,iy,fluids_list(is),istrai)=
     *                   sni(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                      smo(ix,iy,fluids_list(is),istrai)=
     *                   smo(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                    end if !}
                    if(resignore(ix-1,iy-1,1).ne.0) then !{
                      hlp_s(ls)=hlp_s(ls)+
     .                   sni(ix,iy,fluids_list(is),istrai)
                    end if !}
                  end if !}
                end do !}
              end do !}
            end do !}
            do ls=1,nnatmi !{
              hlp_s(ls)=hlp_s(ls)*tflux(istrai)*flux_scale(istrai)
            end do !}

            do j=1,nnlimi+nnstsi !{
              if(j.gt.nnlimi .and. j.ne.l_ddt) then !{
                do is=0,ns-1 !{
                  if(.not.is_neutral(is)) then !{
                    ls=b2espcr(is)
                    if(lstrascl(istrai,ls).gt.0) then !{
                      wldpp(j,fluids_list(is),istrai)=
     =                   wldpp(j,fluids_list(is),istrai)*hlp_scl(ls)
                    end if !}
                  end if !}
                end do !}
                do ls=1,nnatmi !{
                  wldpa(j,ls,istrai)=wldpa(j,ls,istrai)*hlp_scl(ls)
                  hlp_po(ls)=hlp_po(ls)+wldpa(j,ls,istrai)
                end do !}

                do i=1,nnmoli !{
                  ls=lmolscl(i)
                  if(ls.gt.0 .and. ls.le.nnatmi) then !{
                    wldpm(j,i,istrai)=wldpm(j,i,istrai)*hlp_scl(ls)
                  end if !}
                  do ls=1,nnatmi !{
                    hlp_po(ls)=hlp_po(ls)+
     +               wldpm(j,i,istrai)*mlcmp(ls,i)
                  end do !}
                end do !}
              end if !}

              hlp_log=j.le.nnlimi
              if(.not.hlp_log) then !{
                do i=1,ncns !{
                  hlp_log=hlp_log .or. j.eq.lcns(i)+nnlimi
                end do !}
                if(.not.hlp_log) then !{
                  do i=1,ntns !{
                    hlp_log=hlp_log .or. j.eq.ltns(i)+nnlimi
                  end do !}
                end if !}
              end if !}

              do ls=1,nnatmi !{
                if(lstrascl(istrai,ls).gt.0) then !{
                  wldna(j,ls,istrai)=wldna(j,ls,istrai)*hlp_scl(ls)
                  wldra(j,ls,istrai)=wldra(j,ls,istrai)*hlp_scl(ls)
                end if !}
c!###{
c              if(hlp_tt) print '(a,1p,3i4,l3,4e13.5)',
c     ,          ' @@@ istrai,j,ls,log,scl,wldpa,wldna,wldra',
c     ,          istrai,j,ls,hlp_log,hlp_scl(ls),wldpa(j,ls,istrai),
c     ,                        wldna(j,ls,istrai),wldra(j,ls,istrai)
c!###}
                if(hlp_log) then !{
                  hlp_ni(ls)=hlp_ni(ls)+wldna(j,ls,istrai)
                  hlp_no(ls)=hlp_no(ls)+wldra(j,ls,istrai)
                else if(j.eq.l_ddt) then !}{
                  hlp_t(ls)=hlp_t(ls)+wldna(j,ls,istrai)
                end if !}
              end do !}

              do i=1,nnmoli !{
                ls=lmolscl(i)
                if(ls.gt.0 .and. ls.le.nnatmi) then !{
                  wldnm(j,i,istrai)=wldnm(j,i,istrai)*hlp_scl(ls)
                  wldrm(j,i,istrai)=wldrm(j,i,istrai)*hlp_scl(ls)
                end if !}
                do ls=1,nnatmi !{
                  if(hlp_log) then !{
                    hlp_ni(ls)=hlp_ni(ls)+
     +               wldnm(j,i,istrai)*mlcmp(ls,i)
                    hlp_no(ls)=hlp_no(ls)+
     +               wldrm(j,i,istrai)*mlcmp(ls,i)
                  else if(j.eq.l_ddt) then !}{
                    hlp_t(ls)=hlp_t(ls)+
     +               wldnm(j,i,istrai)*mlcmp(ls,i)
                  end if !}
                end do !}
              end do !}
            end do !}
c!###{
            if(dbg_eir_mc.gt.0) then
              do ls=1,nnatmi !{
                s=hlp_po(ls)+hlp_no(ls)-hlp_s(ls)-hlp_ni(ls)
                print '(i6,2i8,3x,1p,8e13.5)',istrai,ls,
     ,                    lstrascl(istrai,ls),tfl(istrai,ls),
     ,                    hlp_pi(ls),hlp_po(ls),hlp_ni(ls),
     ,                    hlp_no(ls),hlp_s(ls),hlp_scl(ls),s
              end do !}
            endif
c!###}

c*** Now rectify the Eirene particle balance for each species -
c*** simply correcting the neutral sources and sinks

            do ls=1,nnatmi !{
c<<<
              if(dbg_eir_mc.gt.0) then
                print *,'@@@ istrai,ls,lstrascl: ',
     ,                            istrai,ls,lstrascl(istrai,ls)
              endif
c>>>
              hlp_scl(ls)=1.0_R8
              if(lstrascl(istrai,ls).eq.ls) then !{

c*** Consider flux of neutrals from plasma recombination on the target
c*** fixed - and adjust other terms to keep the balance

                s=hlp_s(ls)+hlp_t(ls)+hlp_ni(ls)-hlp_no(ls)
                if(s.gt.0.0_R8) hlp_scl(ls)=hlp_po(ls)/s

              else if(lstrascl(istrai,ls).gt.0) then

*** This species in this stratum is born from some interactions of
c*** other species. Try to scale the fluxes to keep the same source
c*** terms in balance.

                s=hlp_po(ls)+hlp_t(ls)+hlp_no(ls)-hlp_ni(ls)
               if(abs(s).gt.0.0_R8 .and.
     .            sign(1.0_R8,s).eq.sign(1.0_R8,hlp_s(ls)))
     >               hlp_scl(ls)=hlp_s(ls)/s
c<<<
               if(dbg_eir_mc.gt.0) then
                 print '(1x,a,1p,5e13.5)','@@@ scl,s,po,no,ni,t: ',
     ,          hlp_scl(ls),s,hlp_po(ls),hlp_no(ls),hlp_ni(ls),hlp_t(ls)
               endif
c>>>

              end if !}
              hlp_scp(ls,istrai)=hlp_scl(ls)
            end do !}

            do ix=1,nx !{
              do iy=1,ny !{
                do is=0,ns-1 !{
                  if(.not.is_neutral(is)) then !{
                    ls=b2espcr(is)
                    if(lstrascl(istrai,ls).gt.0) then !{
                      sni(ix,iy,fluids_list(is),istrai)=
     *                 sni(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                      smo(ix,iy,fluids_list(is),istrai)=
     *                 smo(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                    end if !}
                  end if !}
                end do !}
              end do !}
            end do !}

            do j=1,nnlimi+nnstsi !{
              do ls=1,nnatmi !{
                wldna(j,ls,istrai)=wldna(j,ls,istrai)*hlp_scl(ls)
                wldra(j,ls,istrai)=wldra(j,ls,istrai)*hlp_scl(ls)
                if(lstrascl(istrai,ls).ne.ls) then !{
                  wldpa(j,ls,istrai)=wldpa(j,ls,istrai)*hlp_scl(ls)
                end if !}
              end do !}
              do i=1,nnmoli !{
                ls=lmolscl(i)
                if(ls.gt.0 .and. ls.le.nnatmi) then !{
                  wldnm(j,i,istrai)=wldnm(j,i,istrai)*hlp_scl(ls)
                  wldrm(j,i,istrai)=wldrm(j,i,istrai)*hlp_scl(ls)
                  if(lstrascl(istrai,ls).ne.ls) then !{
                    wldpm(j,i,istrai)=wldpm(j,i,istrai)*hlp_scl(ls)
                  end if !}
                end if !}
              end do !}
            end do !}
          else if(crcstra(istrai).eq.'V' .or.
     .            crcstra(istrai).eq.'C') then !}{

c*** Recombination or gas puffing stratum.
c*** Here only the second pass of correction

            do ls=1,nnatmi !{
              hlp_pi(ls)=0.0_R8
              hlp_po(ls)=0.0_R8
              hlp_ni(ls)=0.0_R8
              hlp_no(ls)=0.0_R8
              hlp_g(ls) =0.0_R8
              hlp_s(ls) =0.0_R8
              hlp_t(ls) =0.0_R8
            end do !}

            do ix=1,nx !{
              do iy=1,ny !{
                if(resignore(ix-1,iy-1,1).ne.0) then !{
                  do is=0,ns-1 !{
                    if(.not.is_neutral(is)) then !{
                      ls=b2espcr(is)
                      hlp_s(ls)=hlp_s(ls)+
     .                 sni(ix,iy,fluids_list(is),istrai)
                    end if !}
                  end do !}
                end if !}
              end do !}
            end do !}
            do ls=1,nnatmi !{
              hlp_s(ls)=hlp_s(ls)*tflux(istrai)*flux_scale(istrai)
            end do !}

            do j=1,nnlimi+nnstsi !{

              hlp_log=j.le.nnlimi
              if(.not.hlp_log) then !{
                do i=1,ncns !{
                  hlp_log=hlp_log .or. j.eq.lcns(i)+nnlimi
                end do !}
                do i=1,ntns !{
                  hlp_log=hlp_log .or. j.eq.ltns(i)+nnlimi
                end do !}
              end if !}

              do ls=1,nnatmi !{
                if(hlp_log) then !{
                  hlp_ni(ls)=hlp_ni(ls)+wldna(j,ls,istrai)
                  hlp_no(ls)=hlp_no(ls)+wldra(j,ls,istrai)
                else if(j.eq.l_ddt) then !}{
                  hlp_t(ls)=hlp_t(ls)+wldna(j,ls,istrai)
                end if !}
              end do !}
              do i=1,nnmoli !{
                ls=lmolscl(i)
                if(ls.gt.0 .and. ls.le.nnatmi) then !{
                  if(hlp_log) then !{
                    hlp_ni(ls)=hlp_ni(ls)+
     +               wldnm(j,i,istrai)*mlcmp(ls,i)
                    hlp_no(ls)=hlp_no(ls)+
     +               wldrm(j,i,istrai)*mlcmp(ls,i)
                  else if(j.eq.l_ddt) then !}{
                    hlp_t(ls)=hlp_t(ls)+
     +               wldnm(j,i,istrai)*mlcmp(ls,i)
                  end if !}
                end if !}
              end do !}
            end do !}

            if(crcstra(istrai).eq.'C') then !{
              do ls=1,nnatmi !{
                kks=ispc_ref(ls)
                if(use_recyceir.gt.0) then
                  recyctmp=recyceir(istrai)
                else
                  recyctmp=recyc(kks,istrai)
                endif
                hlp_g(ls)=tflux(istrai)*gpfc(ls,istrai)*
     *             recyctmp
              end do !}
            end if !}

            do ls=1,nnatmi !{
              hlp_scl(ls)=1.0_R8
              if(lstrascl(istrai,ls).eq.ls) then !{
                s=hlp_s(ls)+hlp_t(ls)
                if(s.eq.0.0_R8) then !{
                  if(hlp_g(ls).gt.0.0_R8) then !{
                    hlp_scl(ls)=(hlp_ni(ls)-hlp_no(ls))/hlp_g(ls)
                  else !}{
                    if(hlp_ni(ls)-hlp_no(ls).ne.0.0_R8)
     -               write(*,'(a,a,a,i3,1x,i3)')
     ,                'eirene_mc: error in scaling - no ',
     ,                'sources and non-zero wall flux. ',
     ,                'stratum, species = ',istrai,ls
                  end if !}
                else !}{
                  s=hlp_g(ls)**2-4.0_R8*(hlp_ni(ls)-hlp_no(ls))*s
                  if(s.gt.0.0_R8) then !{
                    s=sqrt(s)
                    x1=0.5_R8*(hlp_g(ls)+s)/(hlp_s(ls)+hlp_t(ls))
                    x2=0.5_R8*(hlp_g(ls)-s)/(hlp_s(ls)+hlp_t(ls))
                    if(x2.le.0.0_R8) then !{
                      hlp_scl(ls)=x1
                    else if(x1.le.0.0_R8) then !}{
                      hlp_scl(ls)=x2
                    else !}{
                      if(abs(log(x1)).le.abs(log(x2))) then !{
                        hlp_scl(ls)=x1
                      else !}{
                        hlp_scl(ls)=x2
                      end if !}
                    end if !}
                    print '(a,i3,1p,8(1x,e12.5))',
     .                'ank+xpb: ls,s,g,t,ni,no,x1,x2,scl',ls,
     1                 hlp_s(ls),hlp_g(ls),hlp_t(ls),
     .                 hlp_ni(ls),hlp_no(ls),
     2                 x1,x2,hlp_scl(ls)
                  else !}{
                    write(*,'(a,a,i3,1x,i3)')
     ,               'eirene_mc: error in scaling - zero ',
     ,               'discriminant.  stratum, species = ',istrai,ls
                  end if !}
                end if !}
              end if !}
              if(hlp_scl(ls).le.0.0_R8) then !{
                write(*,'(a,a,a,i3,1x,i3,1x,1p,e11.4)')
     ,           'eirene_mc: error in scaling - non-',
     ,           'positive scaling factor. stratum, species, ',
     ,           'factor = ',istrai,ls,hlp_scl(ls)
                hlp_scl(ls)=1.0_R8
              else if(abs(log(hlp_scl(ls))).gt.
     1         abs(log(neut_scl_lim))) then !}{
                write(*,'(a,a,a,i3,1x,i3,1x,1p,e11.4)')
     ,           'eirene_mc: error in scaling - too ',
     ,           'strong variation. stratum, species, ',
     ,           'factor = ',istrai,ls,hlp_scl(ls)
                hlp_scl(ls)=1.0_R8
              end if !}
              hlp_scp(ls,istrai)=hlp_scl(ls)
            end do !}

            do j=1,nnlimi+nnstsi !{
              do ls=1,nnatmi !{
                if(j.ne.l_ddt) then !{
                  wldna(j,ls,istrai)=wldna(j,ls,istrai)/hlp_scl(ls)
                  wldra(j,ls,istrai)=wldra(j,ls,istrai)/hlp_scl(ls)
                else !}{
                  wldna(j,ls,istrai)=wldna(j,ls,istrai)*hlp_scl(ls)
                  wldra(j,ls,istrai)=wldra(j,ls,istrai)*hlp_scl(ls)
                end if !}
              end do !}
              do i=1,nnmoli !{
                ls=lmolscl(i)
                if(ls.gt.0 .and. ls.le.nnatmi) then !{
                  if(j.ne.l_ddt) then !{
                    wldnm(j,i,istrai)=wldnm(j,i,istrai)/hlp_scl(ls)
                    wldrm(j,i,istrai)=wldrm(j,i,istrai)/hlp_scl(ls)
                  else !}{
                    wldnm(j,i,istrai)=wldnm(j,i,istrai)*hlp_scl(ls)
                    wldrm(j,i,istrai)=wldrm(j,i,istrai)*hlp_scl(ls)
                  end if !}
                end if !}
              end do !}
            end do !}

            do ix=1,nx !{
              do iy=1,ny !{
                do is=0,ns-1 !{
                  if(.not.is_neutral(is)) then !{
                    ls=b2espcr(is)
                    if(lstrascl(istrai,ls).gt.0) then !{
                      sni(ix,iy,fluids_list(is),istrai)=
     *                 sni(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                      smo(ix,iy,fluids_list(is),istrai)=
     *                 smo(ix,iy,fluids_list(is),istrai)*hlp_scl(ls)
                    end if !}
                  end if !}
                end do !}
              end do !}
            end do !}
c!###{
            if(dbg_eir_mc.gt.0) then
              print '(i6,2i8,3x,1p,7e13.5)',istrai,ls,
     ,                  lstrascl(istrai,ls),tfl(istrai,ls),
     ,                  hlp_pi(ls),hlp_po(ls),
     ,                  hlp_ni(ls),hlp_no(ls),hlp_s(ls),hlp_scl(ls)
            endif
c!###}
          end if !}
        end do !}

        if(dbg_eir_mc.gt.0 .and. 
     &       main_call .and. mod(ncall,eirene_mod).eq.0) then !{
c*** Print the particle fluxes for each stratum and species
          print '(a)',
     ,     '%%%=============== total fluxes by stratum and '//
     ,     'species after correction ================%%%'
          print '(a25,10(4x,a5,4x))',' Stratum Species lstrascl',
     ,       ' tfl ','wldpi','wldpo','wldni','wldno',
     ,       'Srcpl','puff ',' d/dt','crsrc','error'
          do istrai=1,nstrat !{
            do ls=1,nnatmi !{
              hlp_pi(ls)=0.0_R8
              hlp_po(ls)=0.0_R8
              hlp_ni(ls)=0.0_R8
              hlp_no(ls)=0.0_R8
              hlp_g(ls) =0.0_R8
              hlp_s(ls) =0.0_R8
              hlp_t(ls) =0.0_R8
            end do !}
            do j=1,nnlimi+nnstsi !{
              if(j.gt.nnlimi .and. j.ne.l_ddt) then !{
                do is=0,ns-1 !{
                  if(.not.is_neutral(is)) then !{
                    hlp_pi(b2espcr(is))=hlp_pi(b2espcr(is))+
     +               wldpp(j,fluids_list(is),istrai)
                  end if !}
                end do !}
                do i=1,nnatmi !{
                  hlp_po(i)=hlp_po(i)+wldpa(j,i,istrai)
                  do is=1,nnmoli !{
                    s=mlcmp(i,is)
                    if(s.gt.0.0_R8)
     -               hlp_po(i)=hlp_po(i)+wldpm(j,is,istrai)*s
                  end do !}
                end do !}
              end if !}

              hlp_log=j.le.nnlimi
              if(hlp_log .or. (j.gt.nnlimi .and. j.ne.l_ddt)) then !{
                if(.not.hlp_log) then !{
                  do i=1,ncns !{
                    hlp_log=hlp_log .or. j.eq.lcns(i)+nnlimi
                  end do !}
                  do i=1,ntns !{
                    hlp_log=hlp_log .or. j.eq.ltns(i)+nnlimi
                  end do !}
                end if !}

                if(hlp_log) then !{
                  do i=1,nnatmi !{
                    hlp_ni(i)=hlp_ni(i)+wldna(j,i,istrai)
                    hlp_no(i)=hlp_no(i)+wldra(j,i,istrai)
                    do is=1,nnmoli !{
                      s=mlcmp(i,is)
                      hlp_ni(i)=hlp_ni(i)+wldnm(j,is,istrai)*s
                      hlp_no(i)=hlp_no(i)+wldrm(j,is,istrai)*s
                    end do !}
                  end do !}
                end if !}
              else if(j.eq.l_ddt) then !}{
                do i=1,nnatmi !{
                  hlp_t(i)=hlp_t(i)+wldna(j,i,istrai)
                  do is=1,nnmoli !{
                    s=mlcmp(i,is)
                    hlp_t(i)=hlp_t(i)+wldnm(j,is,istrai)*s
                  end do !}
                end do !}
              end if !}
            end do !}

            do ix=1,nx !{
              do iy=1,ny !{
                if(resignore(ix-1,iy-1,1).ne.0) then !{
                  do is=0,ns-1 !{
                    if(.not.is_neutral(is)) then !{
                      hlp_s(b2espcr(is))=hlp_s(b2espcr(is))+
     +                 sni(ix,iy,fluids_list(is),istrai)
                    end if !}
                  end do !}
                end if !}
              end do !}
            end do !}

            do ls=1,nnatmi !{
              kks=ispc_ref(ls)
              hlp_s(ls)=hlp_s(ls)*tflux(istrai)*flux_scale(istrai)
              t0=hlp_scp(ls,istrai)
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(istrai)
              else
                recyctmp=recyc(kks,istrai)
              endif
              if(crcstra(istrai).eq.'C')
     -           hlp_g(ls)=tflux(istrai)*
     *           gpfc(ls,istrai)*recyctmp
              s=hlp_po(ls)+hlp_g(ls)+hlp_no(ls)-
     -         (hlp_ni(ls)+hlp_s(ls)+hlp_t(ls))
              print '(i6,2i8,3x,1p,10e13.5)',istrai,ls,
     ,         lstrascl(istrai,ls),tfl(istrai,ls),
     ,         hlp_pi(ls),hlp_po(ls),hlp_ni(ls),hlp_no(ls),
     ,         hlp_s(ls),hlp_g(ls),hlp_t(ls),t0,s
              end do !}
          end do !}
          print '(a)',
     ,     '%%%============================================'//
     ,     '=========================================%%%'
        end if !}

c!###{
c        write(hlp_frm1,'(a,i2,a)') '(3x,',nstrat,'i12)'
c        write(hlp_frm2,'(a,i2,a)') '(1x,a5,2x,1p,',nstrat,'e12.4)'
c        do j=1,lns_hlp !{
c          l=lis_hlp(j)
c          print *,'%%%=============== surface ',l
c          write(*,hlp_frm1) (i,i=1,nstrat)
c          write(*,hlp_frm2) 'wldnm',
c     ,                                    (wldnm(l,1,i),i=1,nstrat)
c          write(*,hlp_frm2) 'wldrm',
c     ,                                    (wldrm(l,1,i),i=1,nstrat)
c          write(*,hlp_frm2) 'wldna',
c     ,                                    (wldna(l,1,i),i=1,nstrat)
c          write(*,hlp_frm2) 'wldra',
c     ,                                    (wldra(l,1,i),i=1,nstrat)
c          print *
c        end do !}
c!###}
      endif !ank_mods


csw 26feb2014 accumulate source terms for SOLPS4.3-like averaging procedure
c             (this is new and similar to Kotov's SOLPS4 style, it is NOT IYSs
c              style)
      if(ncall .eq. 0 .and. .not.allocated(sni_av)) then
        allocate(sni_av(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4)))
        allocate(smo_av(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4)))
        allocate(sei_av(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                    size(sei,3)))
        allocate(see_av(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3)))
        sni_av=0.0_R8
        smo_av=0.0_R8
        see_av=0.0_R8
        sei_av=0.0_R8
      endif

      if(neutr_avg .gt. 0) then

       if(main_call.and.mod(ncall,eirene_mod).eq.0) then
        ihist=ihist+1
        if(ihist.eq.neutr_avg+1) then
         sni_av=0.0_R8
         smo_av=0.0_R8
         see_av=0.0_R8
         sei_av=0.0_R8
         dab2_s=0.0_R8
         pfluxa_s=0.0_R8
         pfluxm_s=0.0_R8
         rfluxa_s=0.0_R8
         rfluxm_s=0.0_R8
         ihist = 1
        endif
        sni_av = sni_av + sni
        smo_av = smo_av + smo
        sei_av = sei_av + sei
        see_av = see_av + see
        dab2_s = dab2_s + dab2
        pfluxa_s = pfluxa_s + pfluxa
        pfluxm_s = pfluxm_s + pfluxm
        rfluxa_s = rfluxa_s + rfluxa
        rfluxm_s = rfluxm_s + rfluxm
        dab2 = dab2_s/dble(ihist)
        pfluxa = pfluxa_s/dble(ihist)
        pfluxm = pfluxm_s/dble(ihist)
        rfluxa = rfluxa_s/dble(ihist)
        rfluxm = rfluxm_s/dble(ihist)
       endif
       ih_end=ihist
      else if (nhist.gt.1) then
csw 25feb2014  fix IYS's version of source term averaging (running average)
c fixed SMOOTH SOURCES FROM EIRENE , averaging procedure
       if(main_call.and.mod(ncall,eirene_mod).eq.0) then
        ihist=ihist+1
        if (ihist.eq.nhist) full_history = .true.
        if (ihist.eq.nhist+1) ihist = 1
        write(*,*) 'eirene_mc: EIRENE source averaging ', nhist, ihist
     &      ,full_history
        sni_h(:,:,:,:,ihist) = sni
        smo_h(:,:,:,:,ihist) = smo
        see_h(:,:,:,ihist) = see
        sei_h(:,:,:,ihist) = sei
        do iatm = 1,nnatmi
         pfluxa_h(:,:,iatm,1,ihist) = pfluxa(:,:,iatm,1)
         rfluxa_h(:,:,iatm,1,ihist) = rfluxa(:,:,iatm,1)
         dab2_h(:,:,iatm,1,ihist) = dab2(:,:,iatm,1)
        enddo
        do imol = 1,nnmoli
         pfluxm_h(:,:,imol,1,ihist) = pfluxm(:,:,imol,1)
         rfluxm_h(:,:,imol,1,ihist) = rfluxm(:,:,imol,1)
        enddo
        if (full_history) then
         ih_end = nhist
        else
         ih_end = ihist
        end if
        sni_av = 0.0_R8
        smo_av = 0.0_R8
        see_av = 0.0_R8
        sei_av = 0.0_R8
        dab2_s = 0.0_R8
        pfluxa_s = 0.0_R8      
        pfluxm_s = 0.0_R8      
        rfluxa_s = 0.0_R8      
        rfluxm_s = 0.0_R8      
        do ih = 1, ih_end
         sni_av = sni_av + sni_h(:,:,:,:,ih)
         smo_av = smo_av + smo_h(:,:,:,:,ih)
         see_av = see_av + see_h(:,:,:,ih)
         sei_av = sei_av + sei_h(:,:,:,ih) 
         do iatm = 1,nnatmi
          dab2_s(:,:,iatm,1) = dab2_s(:,:,iatm,1) 
     .                       + dab2_h(:,:,iatm,1,ih)
          pfluxa_s(:,:,iatm,1) = pfluxa_s(:,:,iatm,1) 
     .                         + pfluxa_h(:,:,iatm,1,ih)
          rfluxa_s(:,:,iatm,1) = rfluxa_s(:,:,iatm,1) 
     .                         + rfluxa_h(:,:,iatm,1,ih)
         enddo
         do imol = 1,nnmoli
          pfluxm_s(:,:,imol,1) = pfluxm_s(:,:,imol,1)
     .                         + pfluxm_h(:,:,imol,1,ih)
          rfluxm_s(:,:,imol,1) = rfluxm_s(:,:,imol,1)
     .                         + rfluxm_h(:,:,imol,1,ih)
         enddo
        enddo
        dab2 = dab2_s / dble(ih_end)
        pfluxa = pfluxa_s / dble(ih_end)
        pfluxm = pfluxm_s / dble(ih_end)
        rfluxa = rfluxa_s / dble(ih_end)
        rfluxm = rfluxm_s / dble(ih_end)
       end if
      else if (nrelax.gt.1) then
cxpb Under-relaxation scheme for smoothing of Eirene sources
       if(main_call.and.mod(ncall,eirene_mod).eq.0) then
        ihist=ihist+1
        if (ihist.eq.nrelax) full_history = .true.
        if (ihist.eq.nrelax+1) ihist = 1
        write(*,*) 
     .   'eirene_mc: EIRENE source under-relaxation ', nrelax, ihist
        if (full_history) then
         ih_end = nrelax
        else
         ih_end = ihist
        end if
        relax = rratio(ih_end-1,ih_end)
        sni_s = sni_s*relax + sni/dble(ih_end)
        smo_s = smo_s*relax + smo/dble(ih_end)
        see_s = see_s*relax + see/dble(ih_end)
        sei_s = sei_s*relax + sei/dble(ih_end)
        do iatm = 1,nnatmi
         dab2_s(:,:,iatm,1) = dab2_s(:,:,iatm,1)*relax
     .                      + dab2(:,:,iatm,1)/dble(ih_end)
         pfluxa_s(:,:,iatm,1) = pfluxa_s(:,:,iatm,1)*relax
     .                        + pfluxa(:,:,iatm,1)/dble(ih_end)
         rfluxa_s(:,:,iatm,1) = rfluxa_s(:,:,iatm,1)*relax
     .                        + rfluxa(:,:,iatm,1)/dble(ih_end)
        enddo
        do imol = 1,nnmoli
         pfluxm_s(:,:,imol,1) = pfluxm_s(:,:,imol,1)*relax
     .                        + pfluxm(:,:,imol,1)/dble(ih_end)
         rfluxm_s(:,:,imol,1) = rfluxm_s(:,:,imol,1)*relax
     .                        + rfluxm(:,:,imol,1)/dble(ih_end)
        enddo
        sni_av = sni_s*dble(ih_end)
        smo_av = smo_s*dble(ih_end)
        sei_av = sei_s*dble(ih_end)
        see_av = see_s*dble(ih_end)
        dab2 = dab2_s
        pfluxa = pfluxa_s
        pfluxm = pfluxm_s
        rfluxa = rfluxa_s
        rfluxm = rfluxm_s
       end if
      else
       ih_end=1
       sni_av = sni
       smo_av = smo
       sei_av = sei
       see_av = see
      endif

csw 2012 fixed rescaling procedure
c  otherwise rescaling procedure does not work with internal iterations nstg2>1
c  (sclrtio/sclalph will be reinitialised which should not happen at istg2 > 1)
c  this is part of phase 2 of rescaling process
c 
      if(main_call.and.mod(ncall,eirene_mod).eq.0.and.ank_mods.ne.0)
     .   then

c*** Calculate the basic quantities - i.e. the neutral source strengths
c*** and total plasma sources for each stratum for subsequent rescaling

        do istrai=1, nstrat !{

c*** Accumulators:
c***  sclsrcn : the neutral sources (source strengths)
c***  sclbase : the plasma sources

          do j=0,nnatmi !{
            sclsrcn(j,istrai)=0.0_R8
            sclbase(j,istrai)=0.0_R8
          end do !}
          do is=0,ns-1 !{
            if(.not.is_neutral(is)) then !{
              j=b2espcr(is)
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(istrai)
              else
                recyctmp=recyc(is,istrai)
              endif
              do iy=1,ny !{
                do ix=1,nx !{
                  if(resignore(ix-1,iy-1,1).ne.0) then !{
                    s=sni(ix,iy,fluids_list(is),istrai)*tflux(istrai)*
     *               recyctmp*flux_scale(istrai)
                    if(crcstra(istrai).eq.'V' .and. s.lt.0.0_R8) then !{
                      sclsrcn(j,istrai)=sclsrcn(j,istrai)-s
                      sclsrcn(0,istrai)=sclsrcn(0,istrai)-s
                    else !}{
                      sclbase(j,istrai)=sclbase(j,istrai)+s
                      sclbase(0,istrai)=sclbase(0,istrai)+s
                    end if !}
                  end if !}
                end do !}
              end do !}
            end if !}
          end do !}

          do j=1,nnatmi !{
            is=ispc_ref(j)
            if(use_recyceir.gt.0) then
              recyctmp=recyceir(istrai)
            else
              recyctmp=recyc(is,istrai)
            endif
            if(crcstra(istrai).eq.'E' .or.
     .         crcstra(istrai).eq.'W' .or.
     .         crcstra(istrai).eq.'N' .or.
     .         crcstra(istrai).eq.'S' ) then !{
              sclsrcn(j,istrai)=sclsrcn(j,istrai)+
     +         tfl(istrai,j)*recyctmp
              sclsrcn(0,istrai)=sclsrcn(0,istrai)+
     +         tfl(istrai,j)*recyctmp
            else if(crcstra(istrai).eq.'C'
c     .                          .or. crcstra(istrai).eq.'U'
     .       ) then !}{
              sclsrcn(j,istrai)=sclsrcn(j,istrai)+
     +         tflux(istrai)*gpfc(j,istrai)*recyctmp
            end if !}

            l=lstrascl(istrai,j)
            if(l.gt.nnatmi .or. l.le.0) l=j
            if(sclsrcn(l,istrai).gt.0.0_R8) then !{
              sclbase(j,istrai)=sclbase(j,istrai)/sclsrcn(l,istrai)
            else !}{
              if(sclsrcn(l,istrai).lt.0.0_R8) write(*,*)
     -         '@ neutral source < 0. species,stratum =',l,istrai
              sclbase(j,istrai)=0.0_R8
            end if !}

            sclalph(j,istrai)=1.0_R8
            sclstra(j,istrai)=1.0_R8
          end do !}

          sclsrcn(0,istrai)=0.0_R8
          do j=1,nnatmi !{
            hlp_log=.false.
            do k=1,nnatmi !{
              hlp_log = hlp_log .or. lstrascl(istrai,k).eq.j
            end do !}
            if(hlp_log) sclsrcn(0,istrai)=sclsrcn(0,istrai)+
     +       sclsrcn(j,istrai)
          end do !}
          if(sclsrcn(0,istrai).gt.0.0_R8) then !{
            sclbase(0,istrai)=sclbase(0,istrai)/sclsrcn(0,istrai)
          else !}{
            if(sclsrcn(0,istrai).lt.0.0_R8) write(*,*)
     -       '@ neutral source < 0. species,stratum =',0,istrai
            sclbase(0,istrai)=0.0_R8
          end if !}

          sclalph(0,istrai)=1.0_R8
          sclstra(0,istrai)=1.0_R8
        end do !}
      endif ! ank_mods

      core_flux = 0.0_R8
      if (ionising_core.ne.0) then
        do j = 1, nnlimi+nnstsi
          hlp_cor = .false.
          do i=1,ncns !{
            hlp_cor = hlp_cor .or. j.eq.lcns(i)+nnlimi
          end do !}
          if(hlp_cor) then !{
            do ls=1,nnatmi !{
              core_flux(ls)=core_flux(ls)+(wldna(j,ls,0)-wldra(j,ls,0))
              do is=1,nnmoli !{
                s=mlcmp(ls,is)
                core_flux(ls)=core_flux(ls)+
     .           (wldnm(j,is,0)-wldrm(j,is,0))*s
                if (dbg_eir_mc.ge.12 .and. main_call) then 
                  write(*,'(a,5i4)') 'J, ISTRAI, IATM, IMOL, S',
     .             j, istrai, ls, is, mlcmp(ls,is)
                  write(*,'(a,1p,5(1x,1e14.7))')
     .             'wldna, wldra, wldnm, wldrm, core_flux',
     .              wldna(j,ls,0), wldra(j,ls,0),
     .              wldnm(j,is,0), wldrm(j,is,0), core_flux(ls)
                endif
              end do !}
            end do !}
          end if !}
        end do !}
      endif

csw 2012 now transfer (rescaled) sources into sna0, smo0, she0, shi0,
c   and apply linearisation

      ssni=0.0_R8
      ssmo=0.0_R8
      ssee=0.0_R8
      ssei=0.0_R8
      radarr=0.0_R8
      rad_neutr=0.0_R8
      rdneureg=0.0_R8
      recycireg=0.0_R8
      recycereg=0.0_R8
      sseereg=0.0_R8
      sseireg=0.0_R8
      ssnireg=0.0_R8
      sseicorr=0.0_R8
      src_ni(:,3,:)=0.0_R8
      src_mo(:,3,:)=0.0_R8
      src_ne(:,3)=0.0_R8
      src_ms(:,3)=0.0_R8
      src_ee(:,3)=0.0_R8
      src_ei(:,3)=0.0_R8
      src_es(:,3)=0.0_R8
      src_ne(:,3)=0.0_R8
      smo_tot=0.0_R8
      sni_tot=0.0_R8
      do is = 0, ns-1 !{
        if(.not.is_neutral(is)) then !{
          j=b2espcr(is)
          ifl=fluids_list(is)
          do istrai=1,nstrat !{
            if(use_recyceir.gt.0) then !{
              recyctmp=recyceir(istrai)
            else !}{
              recyctmp=recyc(is,istrai)
            endif !}
            do iy = -1, ny !{
              do ix = -1, nx !{
                smo_tot(ix+1,iy+1,ifl)=smo_tot(ix+1,iy+1,ifl)+
     .           smo_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *           recyctmp*flux_scale(istrai)*sclstra(j,istrai)/
     *           dble(ih_end)
                t00 =
     .           sni_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *           recyctmp*flux_scale(istrai)/dble(ih_end)
                if(crcstra(istrai).ne.'V' .or. t00.gt.0.0_R8) then
                  if (iscltst.eq.0) then
                    t00=t00*sclalph(j,istrai)
                  else
                    t00=t00*sclstra(j,istrai)
                  endif
                endif
                sni_tot(ix+1,iy+1,ifl)=sni_tot(ix+1,iy+1,ifl)+t00
              enddo !}
            enddo !}
          enddo !}
        endif !}
      enddo !}

      do istrai=1,nstrat !{
* by ns
        recyc_hlp=0.0_R8
        do is=0,ns-1 !{
          if(.not.is_neutral(is)) then !{
            ifl=fluids_list(is)
            j=b2espcr(is)
            if(use_recyceir.gt.0) then
              recyctmp=recyceir(istrai)
            else
              recyctmp=recyc(is,istrai)
            endif
            recyc_hlp=max(recyc_hlp,recyctmp)
            do iy=-1,ny !{
              do ix=-1,nx !{
                if(leftix(ix,iy).ne.-2 .and.
     &             rightix(ix,iy).ne.nx+1 .and.
     &             bottomiy(ix,iy).ne.-2 .and.
     &             topiy(ix,iy).ne.ny+1 .and.
     &             region(ix,iy,0).ne.0) then !{
* continuity
                  rf0=RF_NEUT(1)

csw 26feb2014 averaging
c                  t00=sni(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
c     *             recyctmp*flux_scale(istrai)
                  t00=sni_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *             recyctmp*flux_scale(istrai)/
     *             dble(ih_end)

                  if(dpc_fix.lt.0) then
                    write(*,*) 'You are running the code with ' //
     1               'a default option that has been determined ' //
     2               'to be wrong!'
                    write(*,*) 'If you really want to do this, ' //
     1               'you will need to set "eirene_dpc_fix" ' //
     2               'in "b2mn.dat" to "0"'
                    call xertst(dpc_fix.ge.0,
     1               'Invalid setting of "eirene_dpc_fix"')
                    call xerrab ('Error in value of eirene_dpc_fix!')
                  else if(dpc_fix.eq.0) then
                                                  ! original wrong code
                  else if(dpc_fix.eq.1) then
                    t0=t00                        ! dpc fix!
                  else if(dpc_fix.eq.2) then
                    t0=0.0_R8                     ! at least get the wrong code to go wrong safely!
                  endif
                  if(crcstra(istrai).ne.'V' .or. t00.gt.0.0_R8)
     -                                         t0=t00*sclalph(j,istrai)

c                  if (t0.lt.0.0_R8.and.main_call.and.
c     &                crcstra(istrai).ne.'V'.and.output_style.ge.1) then !{
c                    write(6,*) 'neutral particle source < 0! ',t0
c                    write(6,*) 'istrai,is,ix,iy ',istrai,is,ix,iy
c                  else if (t0.gt.0.0_R8.and.main_call.and.
c     &                crcstra(istrai).eq.'V'.and.output_style.ge.1) then !{
c                    write(6,*) 'recombination particle source > 0! ',t0
c                    write(6,*) 'istrai,is,ix,iy ',istrai,is,ix,iy
c                  end if !}

                  if(linearisation.eq.1) then !{
                    if(t0.ge.0.0_R8) then !{
                      sna0(ix,iy,0,is)=sna0(ix,iy,0,is)+t0
                    else !}{
                      sna0(ix,iy,1,is)=sna0(ix,iy,1,is)+t0/na(ix,iy,is)
                    end if !}
                    sna0(ix,iy,0,is)=sna0(ix,iy,0,is)+rf0*abs(t0)
                    sna0(ix,iy,1,is)=sna0(ix,iy,1,is)-
     1               rf0*abs(t0)/na(ix,iy,is)
                  elseif(linearisation.eq.0) then !}{
                    sna0(ix,iy,0,is)=sna0(ix,iy,0,is)+t0
                  end if !}
                  if(crcstra(istrai).ne.'V' .or. t00.gt.0.0_R8)
     >             t0=t00*sclstra(j,istrai)
                  ssni(istrai,is)=ssni(istrai,is)+t0
c*** calculate the integral sources for nregion regions, 3 instants:
c*** before correction (r), after correction (c), after iteration (i)
                  if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                    src_ni(region(ix,iy,0),1,b2espcr(is))=
     .               src_ni(region(ix,iy,0),1,b2espcr(is))+t0
                    src_ni(region(ix,iy,0),2,b2espcr(is))=
     .               src_ni(region(ix,iy,0),2,b2espcr(is))+t0
                    cor_ni(ix,iy,b2espcr(is))=
     .               cor_ni(ix,iy,b2espcr(is))+t0
                    src_ni(0,1,b2espcr(is))=src_ni(0,1,b2espcr(is))+t0
                    src_ni(0,2,b2espcr(is))=src_ni(0,2,b2espcr(is))+t0
                  else
                    src_ni(region(ix,iy,0),3,b2espcr(is))=
     .               src_ni(region(ix,iy,0),3,b2espcr(is))+t0
                    src_ni(0,3,b2espcr(is))=src_ni(0,3,b2espcr(is))+t0
                  endif

csw 21feb2012, add l_neutrad (SOLPS4.3)
                  if(l_neutrad.lt.0) then
                    radarr(ix,iy)=radarr(ix,iy)+t0*rpt(ix,iy,is)*ev
                    rad_neutr=rad_neutr+t0*rpt(ix,iy,is)*ev
                    rdneureg(region(ix,iy,0))=
     .               rdneureg(region(ix,iy,0))+t0*rpt(ix,iy,is)*ev
                    recycereg(region(ix,iy,0))=
     .               recycereg(region(ix,iy,0))+t0*rpt(ix,iy,is)*ev
                  endif
                  ssnireg(b2espcr(is),region(ix,iy,0))=
     .             ssnireg(b2espcr(is),region(ix,iy,0))+t0
* parallel momentum
                  rf0=RF_NEUT(2)
csw 26feb2014 averaging
c                  t0=smo(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
c     *             recyctmp*flux_scale(istrai)*sclstra(j,istrai)
                  t0=smo_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *             recyctmp*flux_scale(istrai)*sclstra(j,istrai)/
     *             dble(ih_end)

                  VTHE=SQRT(TE(IX,IY)/me)
c SCALE VTHE TO GET A VELOCITY CLOSER TO THE ION VELOCITIES
                  VTHE=VTHE/sqrt(am(is)*mp/me)
c
                  if (linearisation.eq.1) then !{
cxpb 2002.01.28 moved to cell-centered velocities
                    smo0(IX,IY,0,IS)=smo0(IX,IY,0,IS)+
     .               t0+RF0*ABS(t0)*UA(IX,IY,IS)/VTHE
                    smo0(IX,IY,1,IS)=smo0(IX,IY,1,IS)-
     .               RF0*ABS(t0)/VTHE
cxpb was                    smo0(IX,IY,0,IS)=smo0(IX,IY,0,IS)+
cxpb was     .               (t0/2.0_R8)+RF0*ABS(t0/2.0_R8)*UA(IX,IY,IS)/VTHE
cxpb was                    smo0(IX,IY,1,IS)=smo0(IX,IY,1,IS)-
cxpb was     .               RF0*ABS((t0/2.0_R8))/VTHE

cxpb was                    smo0(RIGHTIX(IX,IY),RIGHTIY(IX,IY),0,IS)=
cxpb was     1               smo0(RIGHTIX(IX,IY),RIGHTIY(IX,IY),0,IS)+
cxpb was     .               (t0/2.0_R8)+RF0*ABS(t0/2.0_R8)*
cxpb was     1               UA(RIGHTIX(IX,IY),RIGHTIY(IX,IY),IS)/VTHE
cxpb was                    smo0(RIGHTIX(IX,IY),RIGHTIY(IX,IY),1,IS)=
cxpb was     1               smo0(RIGHTIX(IX,IY),RIGHTIY(IX,IY),1,IS)-
cxpb was     .               RF0*ABS((t0/2.0_R8))/VTHE
                  elseif (linearisation.eq.0) then !}{
                    smo0(IX,IY,0,IS)=smo0(IX,IY,0,IS)+t0
                  end if !}
                  SSMO(istrai,IS)=SSMO(istrai,IS)+t0
c*** calculate the integral sources for nregion regions, 3 instants:
c*** before correction (r), after correction (c), after iteration (i)
                  if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                    src_mo(region(ix,iy,0),1,b2espcr(is))=
     .               src_mo(region(ix,iy,0),1,b2espcr(is))+t0
                    src_mo(region(ix,iy,0),2,b2espcr(is))=
     .               src_mo(region(ix,iy,0),2,b2espcr(is))+t0
                    cor_mo(ix,iy,b2espcr(is))=
     .               cor_mo(ix,iy,b2espcr(is))+t0
                    src_mo(0,1,b2espcr(is))=src_mo(0,1,b2espcr(is))+t0
                    src_mo(0,2,b2espcr(is))=src_mo(0,2,b2espcr(is))+t0
                    src_ms(region(ix,iy,0),1)=
     .               src_ms(region(ix,iy,0),1)+t0
                    src_ms(region(ix,iy,0),2)=
     .               src_ms(region(ix,iy,0),2)+t0
                    src_ms(0,1)=src_ms(0,1)+t0
                    src_ms(0,2)=src_ms(0,2)+t0
                  else
                    src_mo(region(ix,iy,0),3,b2espcr(is))=
     .               src_mo(region(ix,iy,0),3,b2espcr(is))+t0
                    src_ms(region(ix,iy,0),3)=
     .               src_ms(region(ix,iy,0),3)+t0
                    src_mo(0,3,b2espcr(is))=src_mo(0,3,b2espcr(is))+t0
                    src_ms(0,3)=src_ms(0,3)+t0
                  endif
                end if !}
              end do !}  ! ix
            end do !}  !iy
          end if !}
        end do !}  ! is
        do iy=-1,ny !{
          do ix=-1,nx !{
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then !{
* electron energy
              rf0=RF_NEUT(3)
csw 26feb2014 averaging
c              t0=see(ix+1,iy+1,istrai)*tflux(istrai)*
c     *         recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)
              t0=see_av(ix+1,iy+1,istrai)*tflux(istrai)*
     *         recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     *         dble(ih_end)

              if (linearisation.eq.1) then !{
                she0(ix,iy,0)=she0(ix,iy,0)+max(0.0_R8,t0)+abs(t0)*rf0
                she0(ix,iy,1)=she0(ix,iy,1)+
     1           (min(0.0_R8,t0)-abs(t0)*rf0)/te(ix,iy)
              elseif (linearisation.eq.0) then !}{
                she0(ix,iy,0)=she0(ix,iy,0)+t0
              end if !}
              ssee(istrai)=ssee(istrai)+t0
c*** calculate the integral sources for nregion regions, 3 instants:
c*** before correction (r), after correction (c), after iteration (i)
              if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                src_ee(region(ix,iy,0),1)=
     .           src_ee(region(ix,iy,0),1)+t0
                src_ee(region(ix,iy,0),2)=
     .           src_ee(region(ix,iy,0),2)+t0
                cor_ee(ix,iy)=cor_ee(ix,iy)+t0
                src_es(region(ix,iy,0),1)=
     .           src_es(region(ix,iy,0),1)+t0
                src_es(region(ix,iy,0),2)=
     .           src_es(region(ix,iy,0),2)+t0
                src_ee(0,1)=src_ee(0,1)+t0
                src_ee(0,2)=src_ee(0,2)+t0
                src_es(0,1)=src_es(0,1)+t0
                src_es(0,2)=src_es(0,2)+t0
              else
                src_ee(region(ix,iy,0),3)=
     .           src_ee(region(ix,iy,0),3)+t0
                src_es(region(ix,iy,0),3)=
     .           src_es(region(ix,iy,0),3)+t0
                src_ee(0,3)=src_ee(0,3)+t0
                src_es(0,3)=src_es(0,3)+t0
              endif
csw 21feb2012, add l_neutrad (SOLPS4.3)
              if(l_neutrad.lt.0) then
                radarr(ix,iy)=radarr(ix,iy)+t0
                rad_neutr=rad_neutr+t0
                rdneureg(region(ix,iy,0))=rdneureg(region(ix,iy,0))+t0
                recycereg(region(ix,iy,0))=recycereg(region(ix,iy,0))
     .                                     +t0
              endif
              sseereg(region(ix,iy,0))=sseereg(region(ix,iy,0))+t0
* ion energy
              rf0=RF_NEUT(4)
csw 26feb2014 averaging
c              t0=sei(ix+1,iy+1,istrai)*tflux(istrai)*
c     1         recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)
              t0=sei_av(ix+1,iy+1,istrai)*tflux(istrai)*
     1         recyc_hlp*flux_scale(istrai)*sclstra(0,istrai)/
     1         dble(ih_end)

              ssei(istrai)=ssei(istrai)+t0
c ? RADARR ?
csw 06mar2012 add l_neutrad
              if(l_neutrad.lt.0) then
                recycireg(region(ix,iy,0))=recycireg(region(ix,iy,0))
     .                                     +t0
              endif
              sseireg(region(ix,iy,0))=sseireg(region(ix,iy,0))+t0
c*** calculate the integral sources for nregion regions, 3 instants:
c*** before correction (r), after correction (c), after iteration (i)
              if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                src_ei(region(ix,iy,0),1)=
     .           src_ei(region(ix,iy,0),1)+t0
                src_ei(region(ix,iy,0),2)=
     .           src_ei(region(ix,iy,0),2)+t0
                cor_ei(ix,iy)=cor_ei(ix,iy)+t0
                src_es(region(ix,iy,0),1)=
     .           src_es(region(ix,iy,0),1)+t0
                src_es(region(ix,iy,0),2)=
     .           src_es(region(ix,iy,0),2)+t0
                src_ei(0,1)=src_ei(0,1)+t0
                src_es(0,1)=src_es(0,1)+t0
                src_ei(0,2)=src_ei(0,2)+t0
                src_es(0,2)=src_es(0,2)+t0
              else
                src_ei(region(ix,iy,0),3)=
     .           src_ei(region(ix,iy,0),3)+t0
                src_es(region(ix,iy,0),3)=
     .           src_es(region(ix,iy,0),3)+t0
                src_ei(0,3)=src_ei(0,3)+t0
                src_es(0,3)=src_es(0,3)+t0
              endif
              t1=0.0_R8
cxpb Correction of Eirene energy sources to B2 heat sources
cxpb Only necessary for SOLPS5.0 equation formulation
              if (BoRiS.eq.0.0_R8.and.no_sseicorr.eq.0) then !{ 
                do is=0,ns-1 !{
                  if(use_recyceir.gt.0) then
                    recyctmp=recyceir(istrai)
                  else
                    recyctmp=recyc(is,istrai)
                  endif
                  if(.not.is_neutral(is)) then !{
                    ifl=fluids_list(is)
cxpb  Moved to cell-center velocities 2002.01.28
cxpb was                 uatmp=0.5_R8*(ua(ix,iy,is)+
cxpb was    1             ua(rightix(ix,iy),rightiy(ix,iy),is))
cxpb  Changed to second-order description to handle very large sources 2007.11.28
cxpb was                 uatmp=ua(ix,iy,is)

csw 06mar2012 include ank_mods scaling into correction
c                    t1=t1-uatmp*smo(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
c     *               recyctmp*flux_scale(istrai)+
c     1               kinrgy(ix,iy,is)*
c     2               sni(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
c     *               recyctmp*flux_scale(istrai)

                    j=b2espcr(is)
                    uatmp = (ua(ix,iy,is)*na(ix,iy,is) + 
     .                       smo_tot(ix+1,iy+1,ifl)*
     .                       dtim*dtmo(is,region(ix,iy,0))/(am(is)*mp))/
     .                      (na(ix,iy,is)+sni_tot(ix+1,iy+1,ifl)*
     .                       dtim*dtco(is,region(ix,iy,0)))
                    dua = uatmp - ua(ix,iy,is)

csw 06mar2012 momentum source
csw 26feb2014 averaging
c                    t2=smo(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
c     *               recyctmp*flux_scale(istrai)
                    t2=smo_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *               recyctmp*flux_scale(istrai)/
     *               dble(ih_end)

                    t2=t2*sclstra(j,istrai)


csw 06mar2012 particle source
csw 26feb2014 averaging
c                    t3=sni(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
c     *               recyctmp*flux_scale(istrai)
                    t3=sni_av(ix+1,iy+1,ifl,istrai)*tflux(istrai)*
     *               recyctmp*flux_scale(istrai)/
     *               dble(ih_end)


csw 26mar2013 (for TESTING)
                    if(iscltst .eq. 0) then
                      if(crcstra(istrai).ne.'V' .or. t3.gt.0.0_R8)
     -                 t3=t3*sclalph(j,istrai)
                    else
                      if(crcstra(istrai).ne.'V' .or. t3.gt.0.0_R8)
     -                 t3=t3*sclstra(j,istrai)
                    endif


csw 06mar2012 total correction
c                    t1=t1-uatmp*t2+kinrgy(ix,iy,is)*t3
cxpb 13feb2015 correction with second-order effects for large sources
                    t1=t1-ua(ix,iy,is)*t2+kinrgy(ix,iy,is)*t3-
     3               (na(ix,iy,is)+t3*dtim*dtco(is,region(ix,iy,0)))/
     *               (dtim*dtei(region(ix,iy,0)))*
     *               0.5_R8*mp*am(is)*dua**2

                  end if !}
                end do !} !is
              end if !}
              sseicorr(istrai)=sseicorr(istrai)+t1
              t0=t0+t1
              if (linearisation.eq.1) then !{
                shi0(ix,iy,0)=shi0(ix,iy,0)+max(0.0_R8,t0)+abs(t0)*rf0
                shi0(ix,iy,1)=shi0(ix,iy,1)+
     1           (min(0.0_R8,t0)-abs(t0)*rf0)/ti(ix,iy)
              elseif (linearisation.eq.0) then !}{
                shi0(ix,iy,0)=shi0(ix,iy,0)+t0
              end if !}
            end if !}
          end do !} ix
        end do !} iy
      end do !} istrai
c
c*** correct radiation losses for the molecule dissociation
c
csw 21feb2012 add test l_neutrad (SOLPS4.3)
      if(l_neutrad.lt.0) then
        do iy=-1,ny
          do ix=-1,nx
            if (leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &          topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &          region(ix,iy,0).ne.0) then
              k=region(ix,iy,0)
              do imol = 1, nnmoli
!pb  eneutrad dimensioned per stratum to allow for short cycling
!pb  only sum over strata used here
!                radarr(ix,iy)=radarr(ix,iy)-edissml(ix+1,iy+1,imol,1)
!                rad_neutr=rad_neutr-edissml(ix+1,iy+1,imol,1)
!                rdneureg(k)=rdneureg(k)-edissml(ix+1,iy+1,imol,1)
!                recycereg(k)=recycereg(k)-edissml(ix+1,iy+1,imol,1)
                radarr(ix,iy)=radarr(ix,iy)-edissml(ix+1,iy+1,imol,0)
                rad_neutr=rad_neutr-edissml(ix+1,iy+1,imol,0)
                rdneureg(k)=rdneureg(k)-edissml(ix+1,iy+1,imol,0)
                recycereg(k)=recycereg(k)-edissml(ix+1,iy+1,imol,0)
              end do
            end if
          end do
        end do
        write (hlp_frm,'(a,i2,a)') '(a,1p,',nnmoli,'(1x,e14.7))'
        if (dbg_eir_mc.ge.2 .and. main_call) write(*,hlp_frm)
     .   'EIRENE_MC. EDISSML : ', 
     .   (SUM(edissml(1:nx,1:ny,imol,0)),imol=1,nnmoli)
      else

csw 21feb2012, calculate neutral radiation directly from EIRENE (SOLPS4.3)
c*** calculate radiation losses directly (at the moment only atoms + edissml)
c
        rad_neutr=0.
        rdneureg=0.
        recycereg=0.
        do iy=-1,ny
          do ix=-1,nx
            if (leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &          topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &          region(ix,iy,0).ne.0) then
              k=region(ix,iy,0)
              do iatm = 1, nnatmi
!pb  18.12.2013
!pb  eneutrad dimensioned per stratum to allow for short cycling
!pb  only sum over strata used here
!                radarr(ix,iy)=radarr(ix,iy)+eneutrad(ix+1,iy+1,iatm,1)
!                rad_neutr=rad_neutr+eneutrad(ix+1,iy+1,iatm,1)
!                rdneureg(k)=rdneureg(k)+eneutrad(ix+1,iy+1,iatm,1)
!                recycereg(k)=recycereg(k)+eneutrad(ix+1,iy+1,iatm,1)
!pb
                 radarr(ix,iy)=eneutrad(ix+1,iy+1,iatm,0)
                 rad_neutr=rad_neutr+eneutrad(ix+1,iy+1,iatm,0)
                 rdneureg(k)=rdneureg(k)+eneutrad(ix+1,iy+1,iatm,0)
                 recycereg(k)=recycereg(k)+eneutrad(ix+1,iy+1,iatm,0)
              enddo
              do imol = 1, nnmoli
!pb  edissml dimensioned per stratum to allow for short cycling
!pb  only sum over strata used here
!                radarr(ix,iy)=radarr(ix,iy)-edissml(ix+1,iy+1,imol,1)
!                rad_neutr=rad_neutr-edissml(ix+1,iy+1,imol,1)
!                rdneureg(k)=rdneureg(k)-edissml(ix+1,iy+1,imol,1)
!                recycereg(k)=recycereg(k)-edissml(ix+1,iy+1,imol,1)
                radarr(ix,iy)=radarr(ix,iy)-edissml(ix+1,iy+1,imol,0)
                rad_neutr=rad_neutr-edissml(ix+1,iy+1,imol,0)
                rdneureg(k)=rdneureg(k)-edissml(ix+1,iy+1,imol,0)
                recycereg(k)=recycereg(k)-edissml(ix+1,iy+1,imol,0)
              end do
            endif
          enddo
        enddo
        write (hlp_frm,'(a,i2,a)') '(a,1p,',nnatmi,'(1x,e14.7))'
        if (dbg_eir_mc.ge.2 .and. main_call) write(*,hlp_frm)
     .   'EIRENE_MC. ENEUTRAD : ', 
     .   (SUM(eneutrad(1:nx,1:ny,iatm,0)),iatm=1,nnatmi)
      endif
      if (dbg_eir_mc.ge.2 .and. main_call) write(*,'(a,1p,2(1x,e14.7))')
     . 'EIRENE_MC. RADARR, RAD_NEUTR : ', 
     .  SUM(radarr(1:nx,1:ny)), rad_neutr

      if(.not.filedata(ntrc_int)) then
        write(*,'(a,1p,1000g12.4)')
     -      'recycecore,recycesol,recycedivin,recycedivout= ',recycereg
        write(*,'(a,1p,1000g12.4)')
     -      'recycicore,recycisol,recycidivin,recycidivout= ',recycireg
      end if

C CORRECTING THE POTENTIAL ENERGY SPENT FOR RECOMBINATION OF ATOMS INTO MOLECULES
      IF(L_NEUTFLUX.GE.0) THEN
        DO I=1,NNSTSI
          k=eirdiag_nds_ind(I)
          if(k.lt.0) cycle
          IS=I+nlim
          wldnep(IS,0)=0.0_R8
          DO J=1,NNMOLI
            wldnep(IS,0)=wldnep(IS,0)+ewldmr_res(j,k)
          END DO
        END DO
      END IF

cxpb  core recycling into ions of the highest ionisation stage
      if (ionising_core.gt.0) then
        rf0=RF_NEUT(1)
        do iy = -1, ny
          do ix = -1, nx
            do j = 1, nnatmi
              if (bottomiy(ix,iy).eq.-2 .and. 
     &            mod(region(ix,iy,0),4).eq.1) then  
                t0=core_flux(j)*
     *           (gs(topix(ix,iy),topiy(ix,iy),1)/core_surface)
                if(linearisation.eq.1) then !{
                  if(t0.gt.0.0_R8) then !{
                    sna0(ix,iy,0,ishigh(j))=sna0(ix,iy,0,ishigh(j))+t0
                  else !}{
                    sna0(ix,iy,1,ishigh(j))=sna0(ix,iy,1,ishigh(j))+
     1               t0/na(ix,iy,ishigh(j))
                  end if !}
                  sna0(ix,iy,0,ishigh(j))=sna0(ix,iy,0,ishigh(j))+
     1             rf0*abs(t0)
                  sna0(ix,iy,1,ishigh(j))=sna0(ix,iy,1,ishigh(j))-
     1             rf0*abs(t0)/na(ix,iy,ishigh(j))
                elseif(linearisation.eq.0) then !}{
                  sna0(ix,iy,0,ishigh(j))=sna0(ix,iy,0,ishigh(j))+t0
                end if !}
              endif
            enddo
          enddo
        enddo
      else if (ionising_core.lt.0) then
        ib=abs(ionising_core)
        call xertst (ib.le.nbc,
     1   'faulty input ionising_core, abs(ionising_core) > NBC')
        do j = 1, nnatmi
          is=ishigh(j)
          if(bccon(is,ib).ne.13) then
            write(*,*) 'eirene_mc: bccon for is,ib= ',is,ib,
     1       ' <> 13, current value = ',bccon(is,ib)
C           call xerrab ('eirene_mc: bccon <> 13')
          else
            conpar(is,ib,3)=core_flux(j)
            if(output_style.gt.2)
     1       write(*,*) 'eirene_mc: conpar(is,ib,3) = ',
     2       conpar(is,ib,3),' for is,ib= ',is,ib
          endif
        enddo
      endif

      cum_volrec=0.0_R8
      do istrai=1, nstrat
        cum_volrec=cum_volrec+VOLSUMN(istrai)
      enddo
      cum_volrec=-cum_volrec/qe

      if((main_call.and.output_style.eq.1).or.output_style.gt.1) then !{
        do istrai=1, nstrat !{
          write (*,*) 'ISTRATUM = ',istrai
          do is=0,ns-1 !{
            write (*,'(a,i3,1p,2(1x,1e14.7))') 'IS,SSNI,SSMO = ',
     1       is,ssni(istrai,is)*qe,ssmo(istrai,is)
          end do !}
          write (*,'(a,1p,4(1x,1e14.7))') 'SSEE,SSEI = ',
     1     ssee(istrai),ssei(istrai)+sseicorr(istrai),
     2     ssei(istrai),sseicorr(istrai)
        end do !}
        write(*,'(a,1p,6(1x,1e14.7))') 
     .  'CUM:SSEE,SSEI[123],SSNI,SSMO = ',
     1   sum(ssee(1:nstrat)),
     2   sum(ssei(1:nstrat))+sum(sseicorr(1:nstrat)),
     3   sum(ssei(1:nstrat)),sum(sseicorr(1:nstrat)),
     2   sum(ssni(1:nstrat,0:ns-1)),sum(ssmo(1:nstrat,0:ns-1))
        if (ionising_core.ne.0) then
          do j = 1, nnatmi 
            write(*,*) 'Recycled',trim(ion_name(j)),
     1       ' ion flux through core boundary', core_flux(j)
          enddo
        endif
        write(*,*) 'Cumulative_volume_recombination ',cum_volrec
      end if !}

      if (main_call.or.ncall.eq.0) ncall=ncall+1

* ..return
      call subend ()
      return

* ..finalise save
      entry eirene_mc_final(nx,ny,ns,BoRiS)
      if(saved_volrec_used) then
        call write_b2mod_neutrals_save(99,'    ')
        call dealloc_b2mod_neutrals_save
      end if
      if(background_new.eq.0) call copy_background(nx, ny, ns, BoRiS)  !srv 14.07.10
      if(background_new.ne.0) call copy_background_(nx, ny, ns, BoRiS) !srv 14.07.10
      call write_f30(nx, ny)
      call write_f31(nx, ny)
csw 26jan2011
      call eirene_extrab25_cleanup
csw 26feb2014
      if(allocated(sni_av)) then
       deallocate(sni_av)
       deallocate(smo_av)
       deallocate(see_av)
       deallocate(sei_av)
      endif
      if(nhist.gt.1) call dealloc_sum_arrays
      if(nrelax.gt.1 .or. nhist.gt.1) call dealloc_hist_arrays

      return

c======================================================================
      entry neutr_scl(nx, ny, ns)
c======================================================================
      if(dbg_eir_mc.gt.0) write(*,*) 'entry neutr_scl'
      if(ank_mods.ne.0) then
c*** Calculate the basic quantities - i.e. the neutral source strengths
c*** and total plasma sources for each stratum for subsequent rescaling

c!###{
        hlp_tt=.false.
c!###}

c*** First, initialise the accumulators.
c***  sclstra: neutral sources
c***  sclrtio: plasma sources

        do istrai=1, nstrat !{
          do j=0,nnatmi !{
            sclstra(j,istrai)=0.0_R8
            sclrtio(j,istrai)=0.0_R8
          end do !}
        end do !}

c*** Evaluate the new volumetric sources cell by cell
c*** Here, the calculation of t0 and s must be consistent with
c*** SNIC and SNIV evaluation above

        do iy=0,ny-1 !{
          do ix=0,nx-1 !{
            if(resignore(ix,iy,1).ne.0) then !{
              do istrai=1, nstrat !{
                do is=0,ns-1 !{
                  if(use_recyceir.gt.0) then
                    recyctmp=recyceir(istrai)
                  else
                    recyctmp=recyc(is,istrai)
                  endif
                  if(.not.is_neutral(is)) then !{
                    j=b2espcr(is)
                    t0=sni(ix+1,iy+1,fluids_list(is),istrai)*recyctmp*
     *               tflux(istrai)*flux_scale(istrai)
c                s=t0+rf0*abs(t0)*(1.0_R8-na(ix,iy,is)/na0(ix,iy,is))
                    if(t0.ge.0.0_R8) then !{
                      s=t0
                    else !}{
                      s=t0*na(ix,iy,is)/na0(ix,iy,is)
                    end if !}
                    if(crcstra(istrai).eq.'V' .and. t0.lt.0.0_R8) then !{
                      sclstra(j,istrai)=sclstra(j,istrai)-s
c                  sclstra(0,istrai)=sclstra(0,istrai)-s
                    else !}{
                      sclrtio(j,istrai)=sclrtio(j,istrai)+
     +                 s*sclalph(j,istrai)
                      sclrtio(0,istrai)=sclrtio(0,istrai)+
     +                 s*sclalph(j,istrai)
                    end if !}
                  end if !}
                end do !}
              end do !}
            end if !}
          end do !}
        end do !}

c*** Calculate the new neutral sources for each stratum. Only the
c*** surface sources should be counted here - the volumetric
c*** recombination is already evaluated above.

        hlp_pff=.true.
        do istrai=1, nstrat !{
          do is=0,ns-1 !{
            if(.not.is_neutral(is)) then !{
              j=b2espcr(is)
              hlp_st=0.0_R8
              if(crcstra(istrai).eq.'W') then !{
                do iy=1,rc_list_size(istrai) !{
                  hlp_st=hlp_st-
     .    fna(rightix(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),
     .        rightiy(rc_list_x(iy,istrai),rc_list_y(iy,istrai)),0,is)
                end do !}
              else if(crcstra(istrai).eq.'E') then !}{
                do iy=1,rc_list_size(istrai) !{
                  hlp_st=hlp_st+
     .             fna(rc_list_x(iy,istrai),rc_list_y(iy,istrai),0,is)
                end do !}
              else if(crcstra(istrai).eq.'N') then !}{
                do ix=1,rc_list_size(istrai) !{
                  hlp_st=hlp_st+fna(rc_list_x(ix,istrai),
     .                              rc_list_y(ix,istrai),1,is)
                end do !}
              else if(crcstra(istrai).eq.'S') then !}{
                do ix=1,rc_list_size(istrai) !{
                  hlp_st=hlp_st-fna(topix(rc_list_x(ix,istrai),
     .                                    rc_list_y(ix,istrai)),
     .                              topiy(rc_list_x(ix,istrai),
     .                                    rc_list_y(ix,istrai)),1,is)
                end do !}
              else if(crcstra(istrai).eq.'C') then !}{
                if(hlp_pff(j,istrai)) then !{
c                if(gpfc(j,istrai).gt.0.0_R8)                    !###
c                  write(*,*) 'istrai,j,is,tflux,gpfc= ',        !###
c                  istrai,j,is,tflux(istrai),gpfc(j,istrai)      !###
                  hlp_st=tflux(istrai)*gpfc(j,istrai)
                  hlp_pff(j,istrai)=.false.
                end if !}
c            else if(crcstra(istrai).eq.'U') then !}{
c              hlp_st=userfluxparm(istrai,1)   ! this is possibly wrong!
              end if !}
              if(use_recyceir.gt.0) then
                recyctmp=recyceir(istrai)
              else
                recyctmp=recyc(is,istrai)
              endif
              sclstra(j,istrai)=sclstra(j,istrai)+
     1         hlp_st*recyctmp
c            sclstra(0,istrai)=sclstra(0,istrai)+hlp_st*recyctmp
            end if !}
          end do !}
        end do !}
c!###{
      if (dbg_eir_mc.gt.0) then
        print *,'@@=============== before correction factors =========='
        print '(a16,20(4x,a7,2x))',
     ,   ' Stratum Species','sclbase','sclstra',
     ,                      'sclrtio','sclalph','sclsrcn'
        do istrai=1, nstrat !{
          do j=0,nnatmi !{
            print '(i6,i8,2x,1p,20e13.5)',istrai,j,sclbase(j,istrai),
     ,         sclstra(j,istrai),sclrtio(j,istrai),sclalph(j,istrai),
     ,                                             sclsrcn(j,istrai)
          end do !}
        end do !}
        print *,'@@===================================================='
      endif
c!###}

c*** Finally, calculate the scaling factors for the plasma sources and
c*** neutral fluxes.
c*** Indeed, and equation F(sclalph)=0 should be solved for each
c*** species and stratum, where F = sclrtio - 1.

        do istrai=1, nstrat !{
          sclstra(0,istrai)=0.0_R8
          do j=1,nnatmi !{
            hlp_stra(j)=sclstra(j,istrai)
c!###{
c          hlp_pr= istrai.eq.7 .and. j.eq.3
c          if(hlp_pr) print *,'%# species,stratum,lstrascl= ',
c     ,                                      j,istrai,lstrascl(istrai,j)
c!###}
            hlp_log=.false.
            do k=1,nnatmi !{
              hlp_log = hlp_log .or. lstrascl(istrai,k).eq.j
            end do !}
            if(hlp_log) sclstra(0,istrai)=sclstra(0,istrai)+
     +       sclstra(j,istrai)
          end do !}
          do j=0,nnatmi !{
c          hlp_pr= istrai.eq.7 .and. j.eq.3 !###
            if(j.gt.0) then !{
              l=lstrascl(istrai,j)
              if(l.le.0 .or. l.gt.nnatmi) l=j
            else !}{
              l=0
            end if !}
            if(sclrtio(j,istrai).gt.0.0_R8 .and. 
     &         sclbase(j,istrai).gt.0.0_R8) then !{
              sclrtio(j,istrai)=sclbase(j,istrai)*
     *         sclstra(l,istrai)/sclrtio(j,istrai)
            else !}{
              if(sclrtio(j,istrai).lt.0.0_R8)
     -         write(0,'(a,a,i3,i3,1p,1e14.7)') 
     ,          '@eirene_mc: negative plasma source. ',
     ,          'istrai,ispecies,source = ',istrai,j,sclrtio(j,istrai)
              sclrtio(j,istrai)=1.0_R8
            end if !}
          end do !}

          do j=0,nnatmi !{
            if(j.gt.0 .and. lstrascl(istrai,j).gt.0 .and.
     .       lstrascl(istrai,j).le.nnatmi) then !{

c*** A primitive iteration scheme - to be modified later if necessary

              l=lstrascl(istrai,j)
              if(l.le.0 .or. l.gt.nnatmi) l=j
              if(sclsrcn(l,istrai).gt.0.0_R8) then !{
                sclstra(j,istrai)=hlp_stra(l)/sclsrcn(l,istrai)
              else !}{
                sclstra(j,istrai)=1.0_R8
              end if !}
              sclstra(j,istrai)=max(1.0_R8/neut_scl_lim,
     1         min(neut_scl_lim,sclstra(j,istrai)))
              if(sclbase(j,istrai).gt.0.0_R8) then !{
                sclalph(j,istrai)=sclalph(j,istrai)*sclrtio(j,istrai)
                sclalph(j,istrai)=max(1.0_R8/neut_scl_lim,
     1           min(neut_scl_lim,sclalph(j,istrai)))
              else !}{
                sclalph(j,istrai)=sclstra(j,istrai)
              end if !}
            else if(j.eq.0 .or. lstrascl(istrai,j).lt.0) then !}{
              if(sclsrcn(j,istrai).gt.0.0_R8) then !{
                sclstra(j,istrai)=sclstra(j,istrai)/sclsrcn(j,istrai)
              else !}{
                sclstra(j,istrai)=1.0_R8
              end if !}
              sclalph(j,istrai)=sclstra(j,istrai)
            else !}{
              sclstra(j,istrai)=1.0_R8
              sclalph(j,istrai)=1.0_R8
            end if !}
          end do !}
        end do !}

c!###{
      if (dbg_eir_mc.gt.0) then
        print *,'@@=============== after correction factors ==========='
        print '(a16,20(4x,a7,2x))',
     ,   ' Stratum Species','sclbase','sclstra',
     ,                      'sclrtio','sclalph','sclsrcn'
        do istrai=1, nstrat !{
          do j=0,nnatmi !{
            print '(i6,i8,2x,1p,20e13.5)',istrai,j,sclbase(j,istrai),
     ,         sclstra(j,istrai),sclrtio(j,istrai),sclalph(j,istrai),
     ,                                             sclsrcn(j,istrai)
          end do !}
        end do !}
        print *,'@@===================================================='
      endif
c!###}

      endif ! ank_mods
      return

      contains

      subroutine alloc_hist_arrays(nhistd)      !IYS 28.03.12 {
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_braeir
#endif
      implicit none
      integer nhistd
#ifndef B25_EIRENE
#include "COUPLE/KOPPLDIM.F"
#endif

! to store the time history of sources from EIRENE
      allocate(sni_h(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4), nhistd))
      allocate(smo_h(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4), nhistd))
      allocate(see_h(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3), nhistd))
      allocate(sei_h(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                    size(sei,3), nhistd))
      allocate(pfluxa_h(0:ndxp,0:ndyp,natm,1,nhistd))
      allocate(pfluxm_h(0:ndxp,0:ndyp,nmol,1,nhistd))
      allocate(rfluxa_h(0:ndxp,0:ndyp,natm,1,nhistd))
      allocate(rfluxm_h(0:ndxp,0:ndyp,nmol,1,nhistd))
      allocate(dab2_h(0:ndxp,0:ndyp,natm,1,nhistd))
      call clear_hist_arrays
      end subroutine alloc_hist_arrays

      subroutine dealloc_hist_arrays
      implicit none

      if (.not.allocated(sni_h)) return
      deallocate(sni_h)
      deallocate(smo_h)
      deallocate(see_h)
      deallocate(sei_h)
      deallocate(pfluxa_h)
      deallocate(pfluxm_h)
      deallocate(rfluxa_h)
      deallocate(rfluxm_h)
      deallocate(dab2_h)
      return
      end subroutine dealloc_hist_arrays

      subroutine clear_hist_arrays                                !csw 24feb2014
      implicit none
      sni_h = 0.0_R8
      smo_h = 0.0_R8
      see_h = 0.0_R8
      sei_h = 0.0_R8
      pfluxa_h = 0.0_R8
      pfluxm_h = 0.0_R8
      rfluxa_h = 0.0_R8
      rfluxm_h = 0.0_R8
      dab2_h = 0.0_R8
      end subroutine clear_hist_arrays                            !IYS 28.03.12 }

      subroutine alloc_sum_arrays
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_braeir
#endif
      implicit none
#ifndef B25_EIRENE
#include "COUPLE/KOPPLDIM.F"
#endif
      allocate(sni_s(0:(size(sni,1)-1), 0:(size(sni,2)-1),
     .                    size(sni,3), size(sni,4)))
      allocate(smo_s(0:(size(smo,1)-1), 0:(size(smo,2)-1),
     .                    size(smo,3), size(smo,4)))
      allocate(see_s(0:(size(see,1)-1), 0:(size(see,2)-1),
     .                    size(see,3)))
      allocate(sei_s(0:(size(sei,1)-1), 0:(size(sei,2)-1),
     .                    size(sei,3)))
      allocate(pfluxa_s(0:ndxp,0:ndyp,natm,1))
      allocate(pfluxm_s(0:ndxp,0:ndyp,nmol,1))
      allocate(rfluxa_s(0:ndxp,0:ndyp,natm,1))
      allocate(rfluxm_s(0:ndxp,0:ndyp,nmol,1))
      allocate(dab2_s(0:ndxp,0:ndyp,natm,1))
      call clear_sum_arrays
      end subroutine alloc_sum_arrays

      subroutine dealloc_sum_arrays
      implicit none
      deallocate(sni_s)
      deallocate(smo_s)
      deallocate(see_s)
      deallocate(sei_s)
      deallocate(pfluxa_s)
      deallocate(pfluxm_s)
      deallocate(rfluxa_s)
      deallocate(rfluxm_s)
      deallocate(dab2_s)
      return
      end subroutine dealloc_sum_arrays

      subroutine clear_sum_arrays                                !csw 24feb2014
      implicit none
      sni_s = 0.0_R8
      smo_s = 0.0_R8
      see_s = 0.0_R8
      sei_s = 0.0_R8
      pfluxa_s = 0.0_R8
      pfluxm_s = 0.0_R8
      rfluxa_s = 0.0_R8
      rfluxm_s = 0.0_R8
      dab2_s = 0.0_R8
      end subroutine clear_sum_arrays                            !IYS 28.03.12 }

      end subroutine eirene_mc
