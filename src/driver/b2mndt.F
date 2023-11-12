*-----------------------------------------------------------------------
*.specification

      subroutine b2mndt (nout, nCv, nFc, nVx, ns,
     . ismain, ismain0, nscx, nscxmax, iscx, itim, dtim, ntim, switch,
     . geo, mpg, st, st_ext, st_avg, ierr)
      use b2mod_types
      use b2mod_constants
      use b2mod_diag
     &   , only : nnatmi, nnmoli, alloc_b2mod_diag
      use b2mod_time
      use b2mod_file
      use b2mod_b2cmpa
      use b2mod_user_namelist
      use b2mod_boundary_namelist
      use b2mod_numerics_namelist
      use b2mod_equation_sources                                        !srv 26.02.18
      use b2mod_boundary_sources                                        !srv 26.02.18
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_running_average
      use b2mod_batch_average
     & , only : batch_av_all
      use b2mod_switches
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_transport_namelist
     & , only : parm_hce, parm_hci, parm_dna, parm_dpa, parm_vla,
     &          parm_vsa, parm_alf, parm_sig
      use b2mod_recycle
     & , only : int0r, int0l, int3r, int3l, int6r, int6l, int2r, int2l,
     &          int5r, int5l, int1r, int1l, int4r, int4l, fna_mol
      use b2mod_b2cmpt
     & , only : cfvla, cflim, cfvsa, cfalf, cfdpa, cfsig, cfdna,
     &          cfhce, cfhci
      use b2mod_input_profile
     & , only : ndata, no_pflux, tr_ip_new_files, region_flags, tdata,
     &          tr_ip_elm_count, no_div, addspec, nss, nrr, nkind_coeff,
     &          nkind_data, poloidal_scaling, transport_ip_filename,
     &          nscale, transport_ip_time_switch, transport_ip_time_mod
      use b2mod_tallies
     & , only : rrahireg, rranareg, rsahireg, rqradreg, rsanareg,
     &          b2bremreg, rqahereg, rramoreg, b2radreg, rqbrmreg,
     &          rsamoreg
      use b2us_feedback
     & , only : fb_target, fb_prev, fb_current, cum_volrec, fb_const,
     &          charge_frac, saved_fb_actuator, fb_rescale, dt_prev,
     &          feedback_namelist_used
      use b2mod_math
     & , only : cutlo, cutll, b2mod_math_initialised, small_r4_constant
#endif
      use b2mod_ad
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nout(0:10), nCv, nFc, nVx, ns, ismain,
     *  ismain0, nscx, nscxmax, iscx(0:nscxmax-1), itim, ntim
      integer, save :: ntim_batch = 500
      integer, save :: ibatch_av_all = 0
      real (kind=R8) ::
     *  dtim
      type (geometry), intent (in) :: geo
      type (mapping), intent (inout) :: mpg
      type (switches), intent (inout) :: switch
*   ..input/output arguments
      type (B2State), intent (inout) :: st
      type (B2StateExt), intent (inout) :: st_ext
      type (B2Average), intent (inout) :: st_avg
*   ..output arguments (unspecified on entry)
      integer ierr
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MNDT performs one implicit timestep for the B2 system of
*     equations.
*
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one timestep for the B2 equations.
*     An outline of the procedure follows.
*     ..test input arguments
*     ..compute auxiliary quantities
*     ..prepare source computation
*     ..do while ()
*      ..compute log-log linearised rate coefficients
*      ..do while ()
*       ..compute source linearisation
*       ..do while ()
*        ..perform one inner iteration
*        ..re-compute auxiliary quantities
*        ..produce monitoring output
*       ..enddo
*      ..enddo
*     ..enddo
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: ramp_slow, ier0, iCv
      integer, save :: ntim_step_out = 1                                 !srv 31.10.07
      integer, save :: use_b2srst = 1                                    !srv 13.07.12
      integer, save :: boundary_sources=0, equation_sources=0            !srv 26.02.18
      save ramp_slow                                                     !xpb
      logical domoqt, wrong_flow                                         !xpb
      logical quit_residual, main_call
      integer istg(0:2), is, k
      real (kind=R8) :: nasum, na0sum
      real (kind=R8) ::
     &  rz(nCv), wrk0(nCv), wrkf(nFc)                                    !iyv 06.02.14
#ifdef WG_TODO
      logical domoit, ramp_call
#endif
      logical solving(4)
!WG_TODO: add some sources locally here, until a module/data-type for WG
!WG_TODO: is available
      real (kind=R8) ::                                                  !xpb
     &  snasm(nCv,0:1,0:ns-1),                                           !xpb
     &  smosm(nCv,0:3,0:ns-1),                                           !xpb
     &  shesm(nCv,0:3), shism(nCv,0:3)!WG_TODO,                          !xpb

*   ..procedures
      real (kind=R8) ::
     *  b2sasum, samax
      external xertst, ipgeti, ipgetr, b2sasum, samax
      external b2xvsg, b2xvps, b2xpne, b2xpni, b2xpnn,
     &  b2spel, b2spcx, b2sral, b2srdt, b2srst, b2mxzr,
     &  b2mxar, b2mxac, b2mwqt, b2mwit, b2tral, b2srsm
*   ..initialisation
      data ramp_slow/0/
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mndt')
*   ..set internal parameters on first call
      if (ncall_b2mndt.eq.0) then
       if (ntim.eq.0) then
         switch%nstg(0) = 1
         switch%nstg(1) = 1
         switch%nstg(2) = 0
       end if
       call ipgeti ('b2news_ramp_slow', ramp_slow)
       if (switch%no_solve.eq.1) switch%nstg(0:2)=1
       test_residual=switch%nstg_areshe.gt.0.0_R8.or.
     &               switch%nstg_areshi.gt.0.0_R8.or.
     &               switch%nstg_aresco.gt.0.0_R8
       if (switch%use_eirene .eq. 1) then
         call alloc_b2mod_diag(nCv,nFc,ns,mpg%nnreg(0),nstrat,nnatmi,
     &                         switch)
       endif
       call ipgeti ('b2mndt_ntim_step_out', ntim_step_out)              !srv 31.01.07 {
       call xertst (1.le.ntim_step_out,
     &   'faulty internal parameter ntim_step_out')                     !srv 31.01.07 }
       call ipgeti ('b2mndt_use_b2srst', use_b2srst)                    !srv 13.07.12
       call ipgeti ('b2mndr_equation_sources', equation_sources)        !srv 26.02.18 {
       if(equation_sources.ne.0) call read_equation_sources(nCv, ns)
       call ipgeti ('b2mndr_boundary_sources', boundary_sources)
       if(boundary_sources.ne.0) call read_boundary_sources(nCv, ns)    !srv 26.02.18 }
       call ipgeti ('b2mndt_av_ntim_batch', ntim_batch)
       call ipgeti ('b2mndt_av_batch_all', ibatch_av_all)

       if (switch%b2mndt_style.eq.2) then
         write(*,*) 'Running in time-dependent mode'
       endif

       if (switch%b2stbc_boundary_namelist.ge.1) then
!WG_TODO         call initialize_b2mod_boundary(ns,mpg,tim)
       else
         call xerrab('b2mndr -- specify BCs using boundary_namelist')
!WG_TODO         call alloc_b2mod_boundary(nx,ny,ns,nbcd,nnisod)
       endif
      else
       if (switch%b2stbc_boundary_namelist.ge.1) then
!WG_TODO         call catchup_b2mod_boundary(ns,mpg,tim)
       endif
      endif
*   ..test nCv, nFc, ns, ismain
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
      call xertst (0.le.ismain0.and.ismain0.lt.ns,
     . 'invalid main neutral species index ismain0')
      call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     . 'invalid main neutral species ismain0; must be neutral')
*   ..test itim, dtim
      call xertst (0.le.itim.and.0.0_R8.lt.dtim,
     &  'faulty argument itim, dtim')
*   ..test sign of bb
      call b2xvsg (nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
*   ..test sign of vol, hx, hy
      call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*   ..test old state
      call b2xvsg (nCv*ns, st%psnc%na, 1, 'na0', '.gt.')
      call b2xvsg (nCv, st%psnc%te, 1, 'te0', '.gt.')
      call b2xvsg (nCv, st%psnc%ti, 1, 'ti0', '.gt.')
      call b2xvsg (nCv, st%psnc%tn, 1, 'tn0', '.gt.')
*   ..test present state
      call b2xvps (nCv, nFc, ns, st%pl, st%dv)
*   ..test facdrift and fac_ExB (not on first call)                      !xpb
      if (ncall_b2mndt.gt.0) then                                        !xpb
        call xertst (0.0_R8.le.minval(st%dv%facdrift) .and.
     &    1.0_R8.ge.maxval(st%dv%facdrift), 'faulty argument facdrift')
        call xertst (0.0_R8.le.minval(st%dv%fac_ExB) .and.
     &    1.0_R8.ge.maxval(st%dv%fac_ExB), 'faulty argument fac_ExB')
      else
c the main initialisation is in b2news --- unfortunately a few routines
c that need an initial value get called before the first call to b2news
c therefore initialise here
cwdk    For unstructured solver, start with simple reading of start-values
        st%dv%fac_ExB  = switch%facExB_start
        st%dv%facdrift = switch%facdrift_start
        st%dv%fac_vis  = switch%facvis_start

#ifdef WG_TODO
        call b2news_facdrift_fac_ExB_init(nx,ny,ramp_slow)
#endif
      endif
*   ..initialise error parameter
      ierr = 0
      quit_residual = .false.

* ..initialisation
*   ..initialise sources to 0
      call sfill (nCv*2*ns, 0.0_R8, snasm, 1)                        !srv 11.09.09 {
      call sfill (nCv*4*ns, 0.0_R8, smosm, 1)
      call sfill (nCv*4, 0.0_R8, shesm, 1)
      call sfill (nCv*4, 0.0_R8, shism, 1)
      call sfill (nCv*2*ns, 0.0_R8, st%sr%snadt, 1)          !! TODO_B2MNDT check !!
      call sfill (nCv*4*ns, 0.0_R8, st%sr%smodt, 1)
      call sfill (nCv*4, 0.0_R8, st%sr%shedt, 1)
      call sfill (nCv*4, 0.0_R8, st%sr%shidt, 1)
      call sfill (nCv*4, 0.0_R8, st%sr%shndt, 1)


*   ..compute ni0
      call b2xpni (nCv, ns, st%psnc%na, st%psnc%ni)
      call b2xpnn (nCv, ns, st%psnc%na, st%psnc%nn)
*   ..compute nasum, na0sum
      nasum = b2sasum(nCv*ns, st%pl%na, 1)
      na0sum = b2sasum(nCv*ns, st%psnc%na, 1)
*   ..rescale densities
      if (switch%density_control.ne.0) then
        st%dv%ne = st%dv%ne*na0sum/nasum
        st%pl%na = st%pl%na*na0sum/nasum
        st%dv%fna = st%dv%fna*na0sum/nasum
        st%dv%kinrgy = st%dv%kinrgy*na0sum/nasum
        write(*,*) 'Densities rescaled by :',na0sum/nasum
      endif
*   ..compute initial ni
      call b2xpni (nCv, ns, st%pl%na, st%dv%ni)
      call b2xpnn (nCv, ns, st%pl%na, st%dv%nn)
* ..main loop
*   ..outer iteration
      istg(0) = 0
    1 if (istg(0).lt.switch%nstg(0) .and. .not.quit_residual) then
*    ..compute log-log electron rate coefficients
        call b2spel (nCv, ns, ev, st%pl%te, st%dv%ne, st%rt)
*    ..compute log-log charge exchange rate coefficients
        do k = 0, nscx-1
         call b2spcx (nCv, ns, ev, am(iscx(k)), st%pl%ti, st%dv%ne,
     &     st%rt%rlcx(1:nCv,0:1,0:ns-1,k))
        enddo
ctmp.dpc
        if (nscx.ge.1) then
          write(*,*) 'DPC: b2mndt: rlcx'
          do k = 0, nscx-1
            do is=0, ns-1
              write(*,'(a,2i4,1p,2g14.7)') 'is0,is,rlcx_0,rlcx_1 ',
     1            k,is,samax(nCv,st%rt%rlcx(1,0,is,k),1),
     2                 samax(nCv,st%rt%rlcx(1,1,is,k),1)
            enddo
          enddo
        else
          write(*,*) 'DPC: b2mndt: no CX reactions'
        endif
ctmp.dpc
*    ..intermediate iteration
        istg(1) = 0
    2   if (istg(1).lt.switch%nstg(1) .and. .not.quit_residual) then
*   ..calculate neoclassical transport only in intermediate loop
!WG_TODO          newneo=.true.

*   ..one-off calculation of transport coefficients
          call b2xpni (nCv, ns, st%pl%na, st%dv%ni)
          call b2xpnn (nCv, ns, st%pl%na, st%dv%nn)
          call b2xpne (nCv, ns, st%rt%rza, st%pl%na, st_ext%ne,
     .      st%dv%ne)
          call b2xpne (nCv, ns, st%rt%rz2, st%pl%na, st_ext%ne2,
     .      st%dv%ne2)
*   ..compute flux limit coefficients                                    !srv 11.10.02 03.06.03 {
          call b2trql (nCv, nFc, ns, switch, geo, mpg, st%pl, st%dv,
     .      st_ext, st%co%chvemx, st%co%chvimx)
          call b2tral (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,   !srv 23.11.10
     1      switch, geo, mpg, st%pl, st%dv, st%rt, st_ext, st%co)
!WG_TODO         newneo=.true.
          call b2tanml (nCv, nFc, ns, switch, geo, mpg, st%co%csig_an,
     &      st%pl%po, st%dv%fchanml)
          call b2tvspa (nCv, nFc, nVx, ns, switch, geo, mpg, st%pl%ua,
     &      st%co%vsaf_cl, st%dv%fac_vis, st%dv%fchvispar)
          call b2tiner (nCv, nFc, nVx, ns, switch, geo, mpg, st%pl%na,
     &      st%pl%ua, st%dv%facdrift, st%dv%fchinert)
#ifdef WG_TODO
          src_ni=0.0_R8
          src_ne=0.0_R8
          src_mo=0.0_R8
          src_ee=0.0_R8
          src_ei=0.0_R8
          src_ms=0.0_R8
          src_es=0.0_R8
          cor_ni=0.0_R8
          cor_ne=0.0_R8
          cor_mo=0.0_R8
          cor_ee=0.0_R8
          cor_ei=0.0_R8
#endif
*     ..compute source linearisation
          wrong_flow = .false.
          main_call = ntim.ne.0
          call b2sral (nCv, nFc, nVx, ns,
     &     nscx, nscxmax, iscx, ismain, ismain0, dtim,
     &     switch, geo, mpg, st, st_ext, st_avg, wrong_flow, main_call)
          if(wrong_flow) write(*,*)
     &         'b2mndt: wrong_flow returned from b2sral'
*     ..smooth out sources
!WG_TODO          call b2srsm (nx, ny, ns,
!WG_TODO     &     dtim, vol, na, ua, te, ti, ne, ni,
!WG_TODO     &     sna, smo, she, shi, sne, .true.)                              !srv 11.09.09 {
*     ..stabilise the source coefficients                                !srv 07.07.21 {
          if(switch%b2mndt_style.eq.2 .and. use_b2srst.eq.1)
     &     call b2srst (nCv, ns, switch, st%pl%na, st%pl%ua,            !srv 01.07.08
     &      st%pl%te, st%pl%ti, st%pl%tn, st%pl%po, st%dv%ne, st%dv%ni,
     &      st%dv%nn, st%pl%kt, st%pl%zt, st%sr)
*       (smq is inherently stable.)                                      !srv 07.07.21 }

*     ..include contributions due to timestep
          if (switch%b2mndt_style.ne.2) then
*     .. contributions due to timestep when style.eq.2 will be done in separate subroutines (b2scdt, b2smdt & b2shdt)
            call b2srdt (nCv, ns, dtim, switch, geo, mpg,
     &       st%psnc%na, st%psnc%ua, st%psnc%te, st%psnc%ti,
     &       st%psnc%tn, st%psnc%ne, st%psnc%ni, st%psnc%nn,
     &       st%psnc%kinrgy,
     &       st%psnc%kt, st%psnc%zt, st%pl%na, st%pl%ua, st%pl%te,
     &       st%pl%ti, st%pl%tn, st%dv%ne, st%dv%ni, st%dv%nn,
     &       st%dv%kinrgy, st%pl%kt, st%pl%zt, st%sr, .false.)
          endif
          if (switch%b2mndt2_backward_comp.eq.1) then !! TODO_B2MNDT is this correct in 3.0.9 ???
*   ..add to previous contributions
            call b2saxpy (nCv*2*ns, 1.0_R8, st%sr%snadt,1, st%sr%sna,1)
            call b2saxpy (nCv*4*ns, 1.0_R8, st%sr%smodt,1, st%sr%smo,1)
            call b2saxpy (nCv*4, 1.0_R8, st%sr%shedt,1, st%sr%she,1)
            call b2saxpy (nCv*4, 1.0_R8, st%sr%shidt,1, st%sr%shi,1)
            call b2saxpy (nCv*4, 1.0_R8, st%sr%shndt,1, st%sr%shn,1)
          endif

          if(equation_sources.ne.0) then                                 !srv 26.02.18 {
            do iCv = 1,nCv
              st%sr%sna(iCv,0:1,0:ns-1) = st%sr%sna(iCv,0:1,0:ns-1)+
     &          geo%cvVol(iCv)*art_sna(iCv,0:1,0:ns-1)
              st%sr%smo(iCv,0:3,0:ns-1) = st%sr%smo(iCv,0:3,0:ns-1)+
     &          geo%cvVol(iCv)*geo%cvHz(iCv)*art_smo(iCv,0:3,0:ns-1)
              st%sr%she(iCv,0:3) = st%sr%she(iCv,0:3)+
     &          geo%cvVol(iCv)*art_she(iCv,0:3)
              st%sr%shi(iCv,0:3) = st%sr%shi(iCv,0:3)+
     &          geo%cvVol(iCv)*art_shi(iCv,0:3)
              st%sr%sch(iCv,0:3) = st%sr%sch(iCv,0:3)+
     &          geo%cvVol(iCv)*art_sch(iCv,0:3)
            enddo
          endif                                                          !srv 26.02.18 }
*     ..stabilise the source coefficients
          if(switch%b2mndt_style.ne.2.and.use_b2srst.eq.1)
     &     call b2srst (nCv, ns, switch, st%pl%na, st%pl%ua,            !srv 01.07.08
     &      st%pl%te, st%pl%ti, st%pl%tn, st%pl%po, st%dv%ne, st%dv%ni,
     &      st%dv%nn, st%pl%kt, st%pl%zt, st%sr)
*       (smq is inherently stable.)
*     ..inner iteration
          istg(2) = 0
    3     if (istg(2).lt.switch%nstg(2) .and. .not.quit_residual) then
*      ..basic iterative step
#ifdef WG_TODO
            ramp_call = ramp_slow.eq.0.or.(istg(0)+istg(1)+istg(2)).eq.0
#endif
            if (switch%b2mndt_style.eq.1) then
              call b2news_(nCv, nFc, nVx, ns, nscx, iscx,
     &         nscxmax, ismain, ismain0, dtim, switch,
     &         geo, mpg, st, st_ext, st_avg, ier0)
            elseif (switch%b2mndt_style.eq.2) then
              call b2news_m(nCv, nFc, nVx, ns, nscx, iscx,
     &         nscxmax, ismain, ismain0, dtim, switch,
     &         geo, mpg, st, st_ext, st_avg, ier0)
            else
              call xerrab('b2mndt_style should be 1 or 2')
            endif
            flush(6)
#ifdef WG_TODO
! ier0 seems not appearing anymore in b2news. Additionally, irreducible loops should be avoided for Adjoint AD
            if(ier0.eq.10) then
             ierr=ier0
             goto 99999
            endif
#endif
            call xertst (ier0.eq.0, 'error return from b2news')
*!!      (should develop other error handling)
            ierr = ier0
*      ..re-compute ne, ne2, ni
            call b2xpne (nCv, ns, st%rt%rza, st%pl%na, st_ext%ne,
     .       st%dv%ne)
            call b2xpne (nCv, ns, st%rt%rz2, st%pl%na, st_ext%ne2,
     .       st%dv%ne2)
            call b2xpni (nCv, ns, st%pl%na, st%dv%ni)
            call b2xpnn (nCv, ns, st%pl%na, st%dv%nn)
*      ..zero out residuals on special regions
*        (this in order to avoid their contribution to the norms)
            call b2mxzr (nCv, ns, mpg, st%dv%resco, st%dv%resmo,
     &       st%dv%resmt, st%dv%reshe, st%dv%reshi, st%dv%reshn,
     &       st%dv%respo, st%dv%reskt, st%dv%reszt)
*      ..compute norms of the residuals
            call b2mxar (nCv, ns, switch,
     &       geo, mpg, st%pl, st%dv, st%diag)
            quit_residual=test_residual .and.
     &       ((st%diag%areshe*switch%nstg_areshe).le.
     &        (switch%nstg_areshe*switch%nstg_areshe)
     &        .or.switch%nstg_areshe.eq.0.0_R8).and.
     &       ((st%diag%areshi*switch%nstg_areshi).le.
     &        (switch%nstg_areshi*switch%nstg_areshi)
     &        .or.switch%nstg_areshi.eq.0.0_R8).and.
     &       ((st%diag%aresco(ismain)*switch%nstg_aresco)
     &        .le.(switch%nstg_aresco*switch%nstg_aresco)
     &        .or.switch%nstg_aresco.eq.0.0_R8)
*      ..compute norms of the corrections
            call b2mxac (nCv, ns, st%dv, st%diag)
*      ..decide what output
            domoqt = 0.le.switch%b2mndt_moqtlv.and.
     &       (1.le.switch%b2mndt_moqtlv.or.istg(0).eq.0).and.
     &       (2.le.switch%b2mndt_moqtlv.or.istg(1).eq.0).and.
     &       (3.le.switch%b2mndt_moqtlv.or.istg(2).eq.0)
#ifdef WG_TODO
            domoit = 0.le.switch%b2mndt_moitlv.and.
     &       (1.le.switch%b2mndt_moitlv.or.istg(0).eq.0).and.
     &       (2.le.switch%b2mndt_moitlv.or.istg(1).eq.0).and.
     &       (3.le.switch%b2mndt_moitlv.or.istg(2).eq.0)
#endif
*      ..maybe write quick trace
            if (domoqt) then
              if (mod(ncall_b2mndt,ntim_step_out) .eq. 0) then             !srv 31.01.07
                call b2mwqt (nout(4), nCv, ns, itim, b2mndt_itcnt,
     &           switch, geo, st%pl, st%dv, st%diag)
                flush(nout(4))
              endif
            endif
#ifdef WG_TODO
*      ..maybe detailed monitor
            if (domoit) then
              call b2mwit (nout(5), nx, ny, ns, itim, b2mndt_itcnt,
     &         resco, resmo, resmt, reshe, reshi, respo,
     &         corpa, corua, corut, corte, corti,
     &         corpo)
              flush(nout(5))
            endif
#endif
*      ..write the tracing data
            if (ank_tracing.ne.0) then
              if (ncall_b2mndt.eq.0.and.istg(2).eq.0) then
                call alloc_b2mod_user(geo,mpg,ns,nlim,nnmoli)
                call b2filei(ns,mpg%nnreg(0),switch)
              end if
              if (mod(itim,ank_tracing).eq.1.or.ank_tracing.eq.1)
     >         call b2filew(istg(2),mpg%nnreg(0))
            endif
*      ..update iteration counter
            b2mndt_itcnt = b2mndt_itcnt+1
*      ..end inner iteration
            istg(2) = istg(2)+1
            goto 3
          else if (ntim.eq.0) then ! compute derived quantities
            if(switch%ion_vlct_restrict .ne. 0 ) then                          !iyv 06.02.14
              call b2xprz (nCv, ns, mp, am, st%pl%na, rz, st_ext)
              call b2xppz (nCv, ns, st%dv%ne, st%pl%na, st%pl%te,
     &          st%pl%ti, st%dv%pz, st_ext)
              wrk0 = sqrt(st%dv%pz/rz)*geo%cvBb(:,0)/geo%cvBb(:,3)
              call intface(nCv,nFc,mpg%fcCv,geo%fcVol,wrk0,st%co%cssb)
            else
              st%co%cssb = 0.0_R8
            endif                                                        !iyv 06
            do is=0,ns-1
*    ..compute partial pressure
              call b2xppb (nCv, st%rt%rza(:,is), st%pl%na(:,is),
     &         st%pl%te, st%pl%ti, st%pl%tn, is, st%dv%pa(:,is))         !srv 01.10.99
*    ..compute particle flux
              call b2tfnb (nCv, nFc, nVx, is, ismain,
     &         switch, geo, mpg, st%pl, st%dv, st%co, st%rt, .true.)
            enddo                                                        !srv 01.10.99
*      ..compute electron drift velocities
            call b2tfed (nCv, nFc, nVx, switch, geo, mpg,
     &       st%dv%facdrift, st%dv%fac_ExB, st%dv%ne, st%pl%te,
     &       st%pl%po, st%dv%vedia, st%dv%veecrb)
*      ..compute other fluxes
            call b2xpen (nCv, nFc, ns, mpg,
     &       st%dv%fna, st%dv%fna_53, st%dv%fne_53,
     &       st%dv%fch, st%dv%fhe, st%dv%fhi, st%dv%fhn, st%rt%rpt,
     &       st%pl%te, st%pl%ti, st%pl%tn, st%pl%po, switch%BoRiS,
     &       st%dv%fhm, st%dv%fhp, st%dv%fhj, st%dv%fht)
*   ..compute ue, electron parallel velocity                             !srv 05.03.11 {
            if (switch%b2sigp_style.eq.2) then
              call b2xpve (nCv, nFc, ns, switch, geo, mpg, qe,
     &          st%rt%rza, st%rt%rz2, st%co%alfx_c, st%co%sigx_c,
     &          st%dv%fch_pi_c, st%pl%po, st%dv%ne,
     &          st%pl%te, st%pl%na, st%pl%ua, st_ext, st%dv%ue)
            else
              wrkf = -st%dv%fch_p(:,0)/(qe*(geo%fcPbs + b2news_cutlo))
              call intcell(nFc,nCv,mpg,mpg%intcellP,wrkf,st%dv%ue)
              do is=0,ns-1
                st%dv%ue = st%dv%ue +
     .                     st%rt%rza(:,is)*st%pl%ua(:,is)*st%pl%na(:,is)
              enddo
              st%dv%ue = st%dv%ue + st_ext%ue*st_ext%ne
              st%dv%ue = st%dv%ue/st%dv%ne                               !srv 05.03.11
            endif
*   ..compute ceqp                                                       !srv 03.07.03 {
            call b2xpnm (nCv, ns, st%rt%rz2, am, st%pl%na, wrk0)
            call b2treq (nCv, switch, geo%cvVol,
     &       st%pl%te, st%pl%ti, st%dv%ne, st%dv%ni,
     &       wrk0, st%dv%lnlam, st_ext, st%co%ceqp)                      !srv 03.07.03 }

*      ..compute energy sources
            solving = .false.
            call b2npht (nCv, nFc, nVx, ns, switch, geo, mpg,
     &       ncall_b2news_, ismain, solving, solveee, solveei,           !srv 22.05.18
     &       solveen, solveet, solvepo, solvemo, solvekt, solvezt,
     &       switch%b2mndt_rxf,
     &       st%pl, st%dv, st%co, st%rt, st%sr, st%srw, st_ext, ier0)
          endif
*     ..end intermediate iteration
          istg(1) = istg(1)+1
          goto 2
        endif
*    ..end outer iteration
        istg(0) = istg(0)+1
        goto 1
      endif

      if (switch%iav_run > 0)
     . call running_average(nCv,ns,st%pl,st%srw,st_avg)
      if (ibatch_av_all > 0) call batch_av_all(itim,ntim_batch,nCv,ns)

* ..return
#ifdef WG_TODO
99999 continue !remove if ierr=10 not supported anymore
#endif
      ncall_b2mndt = ncall_b2mndt+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mndt

      end subroutine b2mndt
