*-----------------------------------------------------------------------
*.specification

      subroutine b2sral (nCv, nFc, nVx, ns,
     &  nscx, nscxmax, iscx, ismain, ismain0, dtim,
     &  switch, geo, mpg, st, st_ext, st_avg, wrong_flow, main_call)
      use b2mod_types
      use b2mod_constants
      use b2mod_tallies
      use b2mod_input_profile
      use b2mod_numerics_namelist
!      use b2mod_boundary_namelist
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_b2cmfs
!      use b2mod_external
      use b2mod_b2_to_astra
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_ad
     & , only : b2sral_elm_count, ncall_b2sral
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_neutrals_namelist
!     & , only : maxw_eff, nstraid, userfluxparm
      use b2mod_recycle
     & , only : int0r, int0l, int3r, int3l, int6r, int6l, int2r, int2l,
     &          int5r, int5l, int1r, int1l, int4r, int4l, fna_mol
      use b2mod_boundary_namelist
     & , only : conpar, enepar, enipar, nbcd, enkpar, potpar, mompar,
     &          lfeedback
      use b2us_feedback
     & , only : fb_target, fb_prev, fb_current, cum_volrec, fb_const,
     &          charge_frac, saved_fb_actuator, fb_rescale, dt_prev,
     &          feedback_namelist_used
      use b2mod_math
     & , only : cutlo, cutll, b2mod_math_initialised, small_r4_constant
#ifndef B25_EIRENE
      use b2mod_eirdiag
     & , only : dab2, lkindi, lkindm, lkindp, pfluxa, rfluxa, tfluxa
#endif
      use b2mod_ad
     & , only : my_out_folder, in_no_of_start_points,
#ifdef WG_DEBUG
     &          my_out_sf,                                              !srv 27.06.21
#endif
     &          ncall_b2stbr_phys, ncall_b2stbr,
     &          ncall_b2stbc_phys, ncall_b2stbc, ncall_b2sqel, nsdmax,
     &          ncall_b2stcx, ncall_b2sqcx, ncall_b2stel, ntstep_b2wall,
     &          ncall_b2tfrn
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns,
     *  nscx, nscxmax, iscx(0:nscxmax-1), ismain, ismain0
      real (kind=R8) :: dtim
      type (switches), intent (inout) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (inout) :: mpg
*   ..input/output arguments
      type (B2State), intent (inout) :: st
      type (B2StateExt), intent (inout) :: st_ext
      type (B2Average), intent (in) :: st_avg
      logical wrong_flow                                                 !xpb
      logical main_call                                                  !dpc
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SRAL provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge.
*     Additionally, an unbalanced electron particle source can be
*     provided.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, k, is, sources_inputfile
      integer, save :: no_fluid_neutral_sources = 0
      character chns*3, chk*1                                           !srv 30.06.08
      real (kind=R8) ::
     *  snatot(0:ns-1),                                                 !srv 18.05.09
!     4  smq0_ion(nCv,0:3,0:ns-1), smq0_rec(nCv,0:3,0:ns-1),
     4  wrk(nCv)!,
!     8  sna0_ion(nCv,0:1,0:ns-1), sna0_rec(nCv,0:1,0:ns-1),
!     9  shi0_ion(nCv,0:3), shi0_rec(nCv,0:3)
      real (kind=R8), allocatable :: sna0_cx(:,:,:), smo0_cx(:,:,:),  ! IYS 06.11.2017
     &  shi0_cx(:,:)
      logical catch_up
      character*32 hlp_frm
      save sources_inputfile
*   ..procedures
      real (kind=R8) ::
     *  smin, smax
      external xertst, sfill, b2saxpy, smin, smax
      external b2xvsg, b2xvps, b2xpfe, b2xpfi, b2sqel, b2sqcx,
     &  b2stbc, b2stbm, b2stbr, b2stel, b2stcx, ipgeti, b2xzdd
*   ..initialisation
      data sources_inputfile/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sral')
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
      call xertst (0.le.ismain0.and.ismain0.lt.ns,
     . 'invalid main neutral species index ismain0')
      call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     . 'invalid main neutral species ismain0; must be neutral')
*   ..inputfile for profile sources
      if (ncall_b2sral.eq.0) then
        call ipgeti ('b2sral_inputfile', sources_inputfile)
        call alloc_input_profile(nCv,ns)
        call ipgeti ('b2sral_no_fluid_neutral_sources',
     .                                       no_fluid_neutral_sources)
      endif
*   ..extensive tests on first few calls
      if (ncall_b2sral.lt.3) then
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
       call b2xvsg (nFc*2, geo%fcVol, 1, 'vol', '.gt.')
*    ..test range of cvQgam, fcQgam
       call xertst (-1.0_R8.lt.smin(nCv,geo%cvQgam(1,1),1).and.
     &   smax(nCv,geo%cvQgam(1,1),1).lt.1.0_R8,
     &   'faulty argument range cvQgam(,1)')
       call xertst (0.0_R8.lt.smin(nCv,geo%cvQgam(1,0),1).and.
     &   smax(nCv,geo%cvQgam(1,0),1).le.1.0_R8,
     &   'faulty argument range cvQgam(,0)')
       call xertst (-1.0_R8.lt.smin(nFc,geo%fcQgam(1,1),1).and.
     &   smax(nFc,geo%fcQgam(1,1),1).lt.1.0_R8,
     &   'faulty argument range fcQgam(,1)')
       call xertst (0.0_R8.lt.smin(nFc,geo%fcQgam(1,0),1).and.
     &   smax(nFc,geo%fcQgam(1,0),1).le.1.0_R8,
     &   'faulty argument range fcQgam(,0)')
*    ..test state
       call b2xvps (nCv, nFc, ns, st%pl, st%dv)
       call b2xvsg (nCv*ns, st%dv%kinrgy, 1, 'kinrgy', '.ge.')
      endif
*   ..compute fne, fni
      call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe,
     & st%rt%rza, st_ext%za, st%dv%fna, st_ext%fa, st%dv%fch, st%dv%fne)
      call b2xpfi (nFc, ns, st%dv%fna, st%dv%fni)

*   ..test facdrift and fac_ExB
      if (ncall_b2sral.eq.0) then                                       !xpb
        if(switch%b2srdt_numerics_namelist.ge.1) then
          call read_b2mod_numerics_namelist(nCv,ns,
     .      switch%nsmin,switch%nsmax,mpg%nnreg,mpg%cvOnClosedSurface,
     .      switch%b2mndt_style,.false.)
          call write_b2mod_numerics_namelist
        endif
        if(numerics_time_mod.gt.0.0_R8) then
          b2sral_elm_count=int(tim/numerics_time_mod)
        else
          b2sral_elm_count=0
        endif
        if(numerics_time_mod.gt.0.0_R8) then
          catch_up=mod(tim,numerics_time_mod).ge.numerics_time_switch
     &      .and.numerics_time_switch.gt.0.0_R8
        else if (numerics_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.numerics_time_switch
        else
          catch_up=.false.
        end if
        do while(catch_up)
          write(*,*) 'b2sral: catching up ...'
          call read_b2mod_numerics_namelist(nCv,ns,
     .      switch%nsmin,switch%nsmax,mpg%nnreg,mpg%cvOnClosedSurface,
     .      switch%b2mndt_style,.false.)
          call write_b2mod_numerics_namelist
          if(numerics_time_mod.gt.0.0_R8) then
            catch_up=mod(tim,numerics_time_mod).ge.numerics_time_switch
     &        .and.numerics_time_switch.gt.0.0_R8
          else if (numerics_time_switch.gt.0.0_R8) then
            catch_up=tim.ge.numerics_time_switch
          else
            catch_up=.false.
          endif
        enddo
      else
        if(numerics_time_mod.gt.0.0_R8) then
          write(*,'(a,1p,3e12.4,2i6)')
     .     'ELM: ',tim,mod(tim,numerics_time_mod),
     1     numerics_time_switch,int(tim/numerics_time_mod),
     2     b2sral_elm_count
          catch_up=
     1    (mod(tim,numerics_time_mod).ge.numerics_time_switch.and.
     1     numerics_time_switch.gt.0.0_R8).or.
     2     int(tim/numerics_time_mod).gt.b2sral_elm_count
        else if (numerics_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.numerics_time_switch
          if(catch_up) write(*,'(a,a)')
     .     'Switching to next numerics namelist file ',
     .      trim(numerics_filename)
        else
          catch_up=.false.
        end if
        if(catch_up) then
          if(numerics_time_mod.gt.0.0_R8)
     >     b2sral_elm_count=int(tim/numerics_time_mod)
          call read_b2mod_numerics_namelist(nCv,ns,
     .     switch%nsmin,switch%nsmax,mpg%nnreg,mpg%cvOnClosedSurface,
     .     switch%b2mndt_style,.false.)
          call write_b2mod_numerics_namelist
        endif
      endif                                                             !xpb
      call xertst (0.0_R8.le.minval(st%dv%facdrift) .and.
     &             1.0_R8.ge.maxval(st%dv%facdrift),
     &            'faulty argument facdrift')
      call xertst (0.0_R8.le.minval(st%dv%fac_ExB) .and.
     &             1.0_R8.ge.maxval(st%dv%fac_ExB),
     &            'faulty argument fac_ExB')

!wdk      if(switch%use_eirene.ne.0 .or. switch%save_f31.ne.0) then     !srv 14.07.10 {
!wdk may want these *_eir-fluxes also without Eirene
*   ..recompute currents
       call b2tanml (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &   st%co%csig_an, st%pl%po, st%dv%ne, st%pl%na, st%rt%rza,
     &   st%dv%fchanml_a, st%dv%fchanml)
       call b2tvspa (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &   st%pl%na, st%pl%ua, st%pl%ti, st%rt%rza, st%co%vsaf_cl,
     &   st%co%vsaf_drho, st%co_ns%vsaf_uAdp_albe,
     &   st%co_ns%vsaf_uBdp_al, st%dv%fac_vis, st%dv%fchvispar_a,
     &   st%dv%fchvispar)
       call b2tiner (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &   st%pl%na, st%pl%ua, st%dv%facdrift,
     &   st%dv%fchinert_a, st%dv%fchinert)
       call b2tfrn (nCv, nFc, nVx, ns, switch, geo, mpg,
     &   st%pl, st%dv, st%co, st%rt, st_ext)
!wdk      endif                                                         !srv 14.07.10 }

* ..compute source coefficients
*   ..compute electron rate coefficients
      call b2sqel (nCv, ns, ismain, switch, ev, st%pl%te, st%rt, st%rtw)
*   ..compute standard contributions due to boundaries
      wrong_flow = .false.                                              !xpb
      if(switch%b2sral_style.eq.0 .or. switch%b2sral_style.eq.2) then   !srv 28.07.08 03.11.17
       if(ncall_b2sral.eq.0) then
         write(*,*) ' style=',switch%b2sral_style,' b2stbc called'
       endif
       call b2stbc (nCv, nFc, nVx, ns, ismain, ismain0, switch,
     &    geo, mpg, st%pl, st%dv, st%co, st%rt, st_ext,
     &    st%srw, st%psnc, st%psnl, wrong_flow, main_call)
*   .. store contributions from BCs
       st%sr%sna = st%srw%sna0
       st%sr%smo = st%srw%smo0
       st%sr%she = st%srw%she0
       st%sr%shi = st%srw%shi0
       st%sr%shn = st%srw%shn0
       st%sr%sch = st%srw%sch0
       st%sr%sne = st%srw%sne0
       st%sr%skt = st%srw%skt0
       st%sr%szt = st%srw%szt0
       if(wrong_flow) write(*,*)
     &  'b2sral [1]: wrong_flow returned from b2stbc'

      else
       call xerrab ('invalid switch b2sral_style')
      endif                                                             !srv 28.07.08 }
#ifdef WG_TODO
      if (switch%use_astra.ne.0) then
        se  =0.0_R8                                                     !srv 28.07.11
        seh =0.0_R8                                                     !srv 11.06.13
        sih =0.0_R8                                                     !srv 11.06.13
        sheNB_astra=0.0_R8                                              !srv 21.06.12 {
        shiNB_astra=0.0_R8
        do iy = -1,ny
          do ix = -1,nx
            if(on_closed_surface(ix,iy)) then                           !srv 23.10.12
              sheNb_astra(ix,iy,0)=vol(ix,iy)*PeNB_astra(iy)*1.0e6_R8
              shiNb_astra(ix,iy,0)=vol(ix,iy)*PiNB_astra(iy)*1.0e6_R8
            endif
          enddo
        enddo
        call b2saxpy (nCv*4, 1.0_R8, sheNb_astra, 1, she, 1)
        call b2saxpy (nCv*4, 1.0_R8, shiNb_astra, 1, shi, 1)             !srv 21.06.12 }
      endif
#endif
*   ..initialise smq to 0
      call sfill (nCv*4*ns, 0.0_R8, st%sr%smq, 1)
*   ..compute special contributions (boundaries and recycling)
      call b2stbr (nCv, nFc, nVx, ns,
     &  nscx, nscxmax, iscx, ismain, dtim,
     &  switch, geo, mpg, st, st_ext, st_avg, main_call)
*   ..zero out contributions in dead cells
      do is = 0, ns-1
        call b2xzdd (nCv, 1, switch, mpg, st%srw%sna0(1,0,is))
        call b2xzdd (nCv, 3, switch, mpg, st%srw%smo0(1,0,is))
      enddo
      call b2xzdd (nCv, 1, switch, mpg, st%srw%sne0)
      call b2xzdd (nCv, 3, switch, mpg, st%srw%she0)
      call b2xzdd (nCv, 3, switch, mpg, st%srw%shi0)
      call b2xzdd (nCv, 3, switch, mpg, st%srw%shn0)
      call b2xzdd (nCv, 3, switch, mpg, st%srw%sch0)
*   ..add to previous contributions
      call b2saxpy (nCv*2*ns, 1.0_R8, st%srw%sna0, 1, st%sr%sna, 1)
      call b2saxpy (nCv*4*ns, 1.0_R8, st%srw%smo0, 1, st%sr%smo, 1)
      call b2saxpy (nCv*2, 1.0_R8, st%srw%sne0, 1, st%sr%sne, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%she0, 1, st%sr%she, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%shi0, 1, st%sr%shi, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%shn0, 1, st%sr%shn, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%sch0, 1, st%sr%sch, 1)
      do is=0,ns-1                                                      !srv 18.05.09 {
        snatot(is)=0.0_R8
        do iCv= 1,nCv
          snatot(is)=snatot(is)+st%sr%sna(iCv,0,is)+
     &     st%sr%sna(iCv,1,is)*st%pl%na(iCv,is)
        enddo
      enddo
      write(hlp_frm,'(a,i3,a)') '(a,1p,',ns,'(1x,1e14.7))'
      if(main_call.or.wrong_flow)                                  !xpb !srv 18.05.09 }
     > write(*,hlp_frm) 'snatot',(snatot(is),is=0,ns-1)
*   ..compute more special contributions
      call b2stbm (nCv, nFc, nVx, ns, ismain, dtim, switch,
     &  geo, mpg, st%pl, st%dv, st_ext, st%srw, main_call)
*   ..zero out contributions in dead cells and core boundary
      do is = 0, ns-1
        call b2xzdd (nCv, 1, switch, mpg, st%srw%sna0(1,0,is))
        call b2xzdd (nCv, 3, switch, mpg, st%srw%smo0(1,0,is))
      enddo
      call b2xzdd (nCv, 1, switch, mpg, st%srw%sne0)
      call b2xzdd (nCv, 3, switch, mpg, st%srw%she0)
      call b2xzdd (nCv, 3, switch, mpg, st%srw%shi0)
      call b2xzdd (nCv, 3, switch, mpg, st%srw%shn0)
      call b2xzdd (nCv, 3, switch, mpg, st%srw%sch0)
*   ..add to previous contributions
      call b2saxpy (nCv*2*ns, 1.0_R8, st%srw%sna0, 1, st%sr%sna, 1)
      call b2saxpy (nCv*4*ns, 1.0_R8, st%srw%smo0, 1, st%sr%smo, 1)
      call b2saxpy (nCv*2, 1.0_R8, st%srw%sne0, 1, st%sr%sne, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%she0, 1, st%sr%she, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%shi0, 1, st%sr%shi, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%shn0, 1, st%sr%shn, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%sch0, 1, st%sr%sch, 1)
*   ..compute sources by inputfiles (heating,...)
      if (sources_inputfile.eq.1) then
        call source_input(mpg,geo,switch,st_ext,ns,
     &   st%srw%sna0,st%srw%smo0,st%srw%she0,st%srw%shi0,st%srw%sch0,
     &   st%srw%sne0,ncall_b2sral)
*   ..zero out contributions in dead cells and core boundary
        do is = 0, ns-1
          call b2xzdd (nCv, 1, switch, mpg, st%srw%sna0(1,0,is))
          call b2xzdd (nCv, 3, switch, mpg, st%srw%smo0(1,0,is))
        enddo
        call b2xzdd (nCv, 1, switch, mpg, st%srw%sne0)
        call b2xzdd (nCv, 3, switch, mpg, st%srw%she0)
        call b2xzdd (nCv, 3, switch, mpg, st%srw%shi0)
        call b2xzdd (nCv, 3, switch, mpg, st%srw%sch0)
*   ..add to previous contributions
        call b2saxpy (nCv*2*ns, 1.0_R8, st%srw%sna0, 1, st%sr%sna, 1)
        call b2saxpy (nCv*4*ns, 1.0_R8, st%srw%smo0, 1, st%sr%smo, 1)
        call b2saxpy (nCv*2, 1.0_R8, st%srw%sne0, 1, st%sr%sne, 1)
        call b2saxpy (nCv*4, 1.0_R8, st%srw%she0, 1, st%sr%she, 1)
        call b2saxpy (nCv*4, 1.0_R8, st%srw%shi0, 1, st%sr%shi, 1)
        call b2saxpy (nCv*4, 1.0_R8, st%srw%sch0, 1, st%sr%sch, 1)
      endif

*   ..compute sources due to electron-atom processes
*     (parallel momentum source must be pure friction)
      call b2stel (nCv, nFc, ns, ismain,
     &  switch, geo, mpg, st%pl, st%dv, st%rt, st%rtw, st%srw) !srv 11.09.09 {
*   ..zero out contributions in dead cells and core boundary
      do is = 0, ns-1
        call b2xzdd (nCv, 1, switch, mpg, st%srw%sna0(1,0,is))
        call b2xzdd (nCv, 3, switch, mpg, st%srw%smq0(1,0,is))
      enddo
      call b2xzdd (nCv, 3, switch, mpg, st%srw%she0)
      call b2xzdd (nCv, 3, switch, mpg, st%srw%shi0)
      call b2xzdd (nCv, 3, switch, mpg, st%srw%shn0)
*   ..add to previous contributions
      call b2saxpy (nCv*2*ns, 1.0_R8, st%srw%sna0, 1, st%sr%sna, 1)
      call b2saxpy (nCv*4*ns, 1.0_R8, st%srw%smq0, 1, st%sr%smq, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%she0, 1, st%sr%she, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%shi0, 1, st%sr%shi, 1)
      call b2saxpy (nCv*4, 1.0_R8, st%srw%shn0, 1, st%sr%shn, 1)
*     (no contribution to sch or sne)
*   ..compute sources due to charge exchange
      rcxnareg=0.0_R8
      rcxhireg=0.0_R8
      rcxmoreg=0.0_R8
      st%srw%rcxna=0.0_R8
      st%srw%rcxhi=0.0_R8
      st%srw%rcxmo=0.0_R8

      if (switch%iout_b2wdat .eq. 4) then  ! IYS 06.11.2017
        if (.not. allocated(shi0_cx)) then
          write(*,*) 'Allocating arrays for CX'
          allocate (shi0_cx(nCv,0:3))
          allocate (sna0_cx(nCv,0:1,0:ns-1))
          allocate (smo0_cx(nCv,0:3,0:ns-1))
          call sfill (nCv*2*ns, 0.0_R8, sna0_cx, 1)
          call sfill (nCv*4*ns, 0.0_R8, smo0_cx, 1)
          call sfill (nCv*4, 0.0_R8, shi0_cx, 1)
        endif
      endif  ! IYS 06.11.2017

      do k = 0, nscx-1
*    ..compute charge exchange rate coefficients
        call b2sqcx (nCv, ns, ev, switch, am(iscx(k)), st%pl%ti,
     &   st%pl%tn, st%rt%rlcx(1:nCv,0:1,0:ns-1,k), st%rtw%rcx)
*    ..compute sources due to charge exchange
*      (parallel momentum source must be pure friction)
        call b2stcx (nCv, nFc, ns, iscx(k), ismain, switch,
     &    geo, mpg, st%pl%na, st%pl%ua, st%pl%ti, st%pl%tn, st%dv%ni,
     &    st%rtw%rcx,
     &    st%srw%sna0, st%srw%smq0, st%srw%shi0, st%srw%shn0,
     &    st%srw%rcxna, st%srw%rcxmo, st%srw%rcxhi)
*   ..zero out contributions in dead cells and core boundary
        do is = 0, ns-1
         call b2xzdd (nCv, 1, switch, mpg, st%srw%sna0(1,0,is))
         call b2xzdd (nCv, 3, switch, mpg, st%srw%smq0(1,0,is))
        enddo
        call b2xzdd (nCv, 3, switch, mpg, st%srw%shi0)
        call b2xzdd (nCv, 3, switch, mpg, st%srw%shn0)
*    ..add to previous contributions
        call b2saxpy (nCv*2*ns, 1.0_R8, st%srw%sna0, 1, st%sr%sna, 1)
        call b2saxpy (nCv*4*ns, 1.0_R8, st%srw%smq0, 1, st%sr%smq, 1)
        call b2saxpy (nCv*4, 1.0_R8, st%srw%shi0, 1, st%sr%shi, 1)
        call b2saxpy (nCv*4, 1.0_R8, st%srw%shn0, 1, st%sr%shn, 1)
*      (no contribution to she, sne or sch)
        if (switch%iout_b2wdat .eq. 4) then  ! IYS 06.11.2017
          write(*,*) 'Summing cx sources over hydrogenic neutrals'
          call b2saxpy (nCv*2*ns, 1.0_R8, st%srw%sna0, 1, sna0_cx, 1)
          call b2saxpy (nCv*4*ns, 1.0_R8, st%srw%smq0, 1, smo0_cx, 1)
          call b2saxpy (nCv*4, 1.0_R8, st%srw%shi0, 1, shi0_cx, 1)
        endif                         ! IYS 06.11.2017
      enddo

      if (switch%iout_b2wdat .eq. 4) then  ! IYS 06.11.2017
        write(*,*) 'Writing CX sources'
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          wrk(:) = sna0_cx(:,0,is)+sna0_cx(:,1,is)*st%pl%na(:,is)
          call my_out_us(70,nCv,0,wrk,'b2stcx_sna_'//chns)
          wrk(:) = shi0_cx(:,0)+shi0_cx(:,1)*st%pl%ti+
     &             shi0_cx(:,2)*st%dv%ni(:,0)+
     &             shi0_cx(:,3)*st%pl%ti*st%dv%ni(:,0)
          call my_out_us(70,nCv,0,wrk,'b2stcx_shi')
          wrk(:) = smo0_cx(:,0,is)+
     &             smo0_cx(:,1,is)*st%pl%ua(:,is)+
     &             smo0_cx(:,2,is)*st%pl%na(:,is)*mp*am(is)+
     &             smo0_cx(:,3,is)*st%pl%ua(:,is)
     &                            *st%pl%na(:,is)*mp*am(is)
          call my_out_us(70,nCv,0,wrk,'b2stcx_smq'//chns)
        enddo
        if (allocated(shi0_cx)) then
          write(*,*) 'Deallocating arrays for cx'
          deallocate (shi0_cx)
          deallocate (sna0_cx)
          deallocate (smo0_cx)
        endif
      endif

      if(switch%b2sral_iout.ne.0) then                                   !srv 30.06.08 {
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          do k=0,1
            write (chk,'(i1)') k
            call my_out_us(70,nCv,0,st%sr%sna(1,k,is),
     &       'b2sral_sna'//chk//chns)
          enddo
          do k=0,3
            write (chk,'(i1)') k
            call my_out_us(70,nCv,0,st%sr%smo(1,k,is),
     &       'b2sral_smo'//chk//chns)
            call my_out_us(70,nCv,0,st%sr%smq(1,k,is),
     &       'b2sral_smq'//chk//chns)
          enddo
        enddo
        do k=0,3
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,st%sr%she(1,k),'b2sral_she'//chk)
          call my_out_us(70,nCv,0,st%sr%shi(1,k),'b2sral_shi'//chk)
          call my_out_us(70,nCv,0,st%sr%shn(1,k),'b2sral_shn'//chk)
          call my_out_us(70,nCv,0,st%sr%sch(1,k),'b2sral_sch'//chk)
#ifdef WG_TODO
          if (switch%use_astra.ne.0) then
            call my_out_us(70,nCv,0,sheNb_astra(1,k),
     &       'b2sral_sheNb_astra'//chk)
            call my_out_us(70,nCv,0,shiNb_astra(1,k),
     &       'b2sral_shiNb_astra'//chk)
          endif
#endif
        enddo
        do k=0,1
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,st%sr%sne(1,k),'b2sral_sne'//chk)
        enddo
      endif                                                              !srv 30.06.08 }

csw 17oct2012, erase sna, smo, smq, (she, shi, sch) for fluid neutrals
      if(no_fluid_neutral_sources .eq. 1) then
       do is=0,ns-1
        if(is_neutral(is)) then
          st%sr%sna(:,:,is) =0.d0
          st%sr%smo(:,:,is) =0.d0
          st%sr%smq(:,:,is) =0.d0
        endif
       enddo
      endif

* ..return
      if(ncall_b2sral.eq.0 .or. main_call) ncall_b2sral = ncall_b2sral+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sral

      end subroutine b2sral
