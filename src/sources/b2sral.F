*-----------------------------------------------------------------------
*.specification

      subroutine b2sral (nCv, nFc, nVx, ns, nxtl, nxtr,
     &  nscx, nscxmax, nsmin, nsmax, iscx, ismain, ismain0, dtim, BoRiS,
     &  geo, mpg, pl, dv, co, rt, rtw, st_ext, sr, srw,
     &  wrong_flow, main_call)
      use b2mod_types
      use b2mod_constants
      use b2mod_tallies
!      use b2mod_sources
!     &     , only: rcxna, rcxhi, rcxmo
      use b2mod_input_profile
      use b2mod_numerics_namelist
!      use b2mod_boundary_namelist
      use b2mod_time
!      use b2mod_anomalous_transport
      use b2mod_b2cmpa
      use b2mod_b2cmfs
!      use b2mod_external
      use b2mod_b2_to_astra
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nxtl, nxtr,
     *  nscx, nscxmax, iscx(0:nscxmax-1), ismain, ismain0, nsmin, nsmax
      real (kind=R8) ::                                                  !xpb
     *  dtim, BoRiS
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
*   ..input/output arguments
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (inout) :: rt
      type (B2RatesWork), intent (inout) :: rtw
      type (B2Coeff), intent (inout) :: co
      type (B2StateExt), intent (inout) :: st_ext
      type (B2Source), intent (inout) :: sr
      type (B2SourceWork), intent (inout) :: srw
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
      integer ncall, iCv, k, is, sources_inputfile
      integer numerics_namelist
      integer use_astra
      integer, save :: elm_count
      integer, save :: style = 2                                        !srv 28.06.08
      integer, save :: iout = 0, iout_b2wdat = 0                        !srv 30.06.08, IYS 06.11.2017
      integer, save :: use_eirene =0, save_f30 = 0, save_f31 = 0        !srv 14.07.10
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
      save ncall, sources_inputfile, numerics_namelist
      save use_astra
*   ..procedures
      real (kind=R8) ::
     *  smin, smax
      external subini, subend, xertst, sfill, b2saxpy, smin, smax
      external b2xvsg, b2xvff, b2xvps, b2xpfe, b2xpfi, b2sqel, b2sqcx,
     &  b2stbc, b2stbm, b2stbr, b2stel, b2stcx, ipgeti, b2xzdd
*   ..initialisation
      data ncall/0/, sources_inputfile/0/, use_astra/0/,
     *  numerics_namelist/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sral')
*   ..test nCv, nFc, ns
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
*   ..test facdrift (not on first call)
      if (ncall.gt.0) then
        call xertst (0.0_R8.le.minval(dv%facdrift) .and.
     &       1.0_R8.ge.maxval(dv%facdrift), 'faulty argument facdrift')
        call xertst (0.0_R8.le.minval(dv%fac_ExB) .and.
     &       1.0_R8.ge.maxval(dv%fac_ExB), 'faulty argument fac_ExB')
       endif
*   ..inputfile for profile sources
      if (ncall.eq.0) then
        call ipgeti ('b2sral_inputfile', sources_inputfile)
        call alloc_input_profile(nCv,ns)
        call ipgeti ('b2sral_style', style)                              !srv 28.07.08
        call ipgeti ('b2sral_iout', iout)                                !srv 30.06.08
        call ipgeti ('b2mndr_astra', use_astra)
        call ipgeti ('b2mndr_eirene', use_eirene)                        !srv 14.07.10 {
        call ipgeti ('eirene_savef30', save_f30)
        call ipgeti ('eirene_savef31', save_f31)                         !srv 14.07.10 }
        call ipgeti ('b2sral_no_fluid_neutral_sources',
     .                                       no_fluid_neutral_sources)
        call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 06.11.17
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
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
       call b2xvps (nCv, nFc, ns, pl, dv)
       call b2xvsg (nCv*ns, dv%kinrgy, 1, 'kinrgy', '.ge.')
      endif
*   ..compute fne, fni
      call b2xpfe (nCv, nFc, ns, st_ext%ns, geo, mpg, qe, 
     &  rt%rza, st_ext%za, dv%fna, st_ext%fa, dv%fch, dv%fne)
      call b2xpfi (nFc, ns, dv%fna, dv%fni)

*   ..test facdrift and fac_ExB
      if (ncall.eq.0) then                                              !xpb
        call ipgeti ('b2srdt_numerics_namelist', numerics_namelist)
        call alloc_b2mod_numerics(nCv,ns,mpg%nnreg)
        if(numerics_namelist.ge.1) then
          call read_b2mod_numerics_namelist(nCv,nsmin,nsmax,mpg%nnreg,
     .      mpg%cvOnClosedSurface)
          call write_b2mod_numerics_namelist
        endif
        if(numerics_time_mod.gt.0.0_R8) then
          elm_count=tim/numerics_time_mod
        else
          elm_count=0
        endif
        catch_up=numerics_time_mod.gt.0.0_R8
        do while(catch_up)
          if(mod(tim,numerics_time_mod).gt.numerics_time_switch)
     1     then
            write(*,*) 'b2srsm: catching up ...'
            call read_b2mod_numerics_namelist(nCv,nsmin,nsmax,mpg%nnreg,
     .        mpg%cvOnClosedSurface)
            call write_b2mod_numerics_namelist
            catch_up=numerics_time_mod.gt.0.0_R8
          else
            catch_up=.false.
          endif
        enddo
      else
        if(numerics_time_mod.gt.0.0_R8) then
          write(*,'(a,1p,3e12.4,2i6)')
     .     'ELM: ',tim,mod(tim,numerics_time_mod),
     1     numerics_time_switch,int(tim/numerics_time_mod),elm_count
          if(mod(tim,numerics_time_mod).gt.numerics_time_switch.or.
     1     int(tim/numerics_time_mod).gt.elm_count) then
            elm_count=tim/numerics_time_mod
            call read_b2mod_numerics_namelist(nCv,nsmin,nsmax,mpg%nnreg,
     .        mpg%cvOnClosedSurface)
            call write_b2mod_numerics_namelist
          endif
        endif
      endif                                                             !xpb
      call xertst (0.0_R8.le.minval(dv%facdrift) .and.
     &     1.0_R8.ge.maxval(dv%facdrift), 'faulty argument facdrift')
      call xertst (0.0_R8.le.minval(dv%fac_ExB) .and.
     &     1.0_R8.ge.maxval(dv%fac_ExB), 'faulty argument fac_ExB')

      if(use_eirene.ne.0 .or. save_f31.ne.0) then                       !srv 14.07.10 {
*   ..recompute currents
       call b2tanml (nCv, nFc, ns, geo, mpg, co%csig_an, pl%po,
     &   dv%fchanml)
       call b2tvspa (nCv, nFc, nVx, ns, geo, mpg, pl%ua, co%vsaf_cl, 
     &   dv%fac_vis, dv%fchvispar)
       call b2tiner (nCv, nFc, nVx, ns, geo, mpg, pl%na, pl%ua,
     &      dv%facdrift, dv%fchinert)
       call b2tfrn (nCv, nFc, nVx, ns, geo, mpg, pl, dv, co, rt, st_ext)
      endif                                                             !srv 14.07.10 }

* ..compute source coefficients
*   ..compute electron rate coefficients
      call b2sqel (nCv, ns, ismain, ev, pl%te, rt, rtw)
*   ..compute standard contributions due to boundaries
      wrong_flow = .false.                                              !xpb
      if(style.eq.0 .or. style.eq.2) then                               !srv 28.07.08 03.11.17
       if(ncall.eq.0) write(*,*) ' style=',style,' b2stbc called'
       call b2stbc (nCv, nFc, nVx, ns, ismain, ismain0, BoRiS, geo,
     &   mpg, pl, dv, co, rt, rtw, st_ext, srw, wrong_flow, main_call)
*   .. store contributions from BCs
       sr%sna = srw%sna0
       sr%smo = srw%smo0
       sr%she = srw%she0
       sr%shi = srw%shi0
       sr%shn = srw%shn0
       sr%sch = srw%sch0
       sr%sne = srw%sne0
       if(wrong_flow) write(*,*)
     &  'b2sral [1]: wrong_flow returned from b2stbc'


#ifdef WG_TODO
      elseif(style.eq.1) then
!! FIXME : Need to add sne in argument list here ???
       if(ncall.eq.0) write(*,*) ' style=',style,' b2stbc_spb called'
       call b2stbc_spb (nx, ny, ns, facdrift, vol, gs, hx, hy, hz, qz,  !srv 09.01.01
     &  pbs, pbshz, bb, na, ua, uadia, vaecrb, vadia, te, ti, po, ne,   !srv 14.10.99 12.05.00 09.01.01 28.10.02
     1  ni, floe_noc, floi_noc, fna, fna_fcor, fnaPSch, fhe, fhi,       !srv 13.10.06
     &  fhePSch, fhiPSch, fne, fni, fch,                                !srv 25.02.00 21.10.02
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,  !srv 21.10.02
     &  sna, smo, she, shi, sch, wrong_flow)                            !srv 09.07.99 01.07.05
       if(wrong_flow) write(*,*)
     &  'b2sral [2]: wrong_flow returned from b2stbc_spb'
#endif
      else
       call xertst (0.le.style.and.style.le.2,                          !srv 03.11.17
     .  'invalid switch b2sral_style')
      endif                                                             !srv 28.07.08 }
#ifdef WG_TODO
      if (use_astra.ne.0) then
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
      call sfill (nCv*4*ns, 0.0_R8, sr%smq, 1)
*   ..compute special contributions (boundaries and recycling)
      call b2stbr (nCv, nFc, nVx, ns, nxtl, nxtr, 
     &  nscx, nscxmax, iscx, ismain, dtim, BoRiS,
     &  geo, mpg, pl, dv, co, rt, rtw, st_ext, srw, main_call)
*   ..zero out contributions in dead cells
      do is = 0, ns-1
        call b2xzdd (nCv, 1, mpg, srw%sna0(1,0,is))
        call b2xzdd (nCv, 3, mpg, srw%smo0(1,0,is))
      enddo
      call b2xzdd (nCv, 1, mpg, srw%sne0)
      call b2xzdd (nCv, 3, mpg, srw%she0)
      call b2xzdd (nCv, 3, mpg, srw%shi0)
      call b2xzdd (nCv, 3, mpg, srw%shn0)
      call b2xzdd (nCv, 3, mpg, srw%sch0)
*   ..add to previous contributions
      call b2saxpy (nCv*2*ns, 1.0_R8, srw%sna0, 1, sr%sna, 1)
      call b2saxpy (nCv*4*ns, 1.0_R8, srw%smo0, 1, sr%smo, 1)
      call b2saxpy (nCv*2, 1.0_R8, srw%sne0, 1, sr%sne, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%she0, 1, sr%she, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%shi0, 1, sr%shi, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%shn0, 1, sr%shn, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%sch0, 1, sr%sch, 1)
      do is=0,ns-1                                                      !srv 18.05.09 {
        snatot(is)=0.0_R8
        do iCv= 1,nCv
          snatot(is)=snatot(is)+sr%sna(iCv,0,is)+
     &     sr%sna(iCv,1,is)*pl%na(iCv,is)
        enddo
      enddo
      write(hlp_frm,'(a,i3,a)') '(a,1p,',ns,'(1x,1e14.7))'
      if(main_call.or.wrong_flow)                                  !xpb !srv 18.05.09 }
     > write(*,hlp_frm) 'snatot',(snatot(is),is=0,ns-1)


*   ..compute more special contributions
      call b2stbm (nCv, nFc, nVx, ns, ismain, dtim,
     &  geo, mpg, pl, dv, co, rt, rtw, st_ext, sr, srw,
     &  main_call)
#ifdef WG_TODO
*   ..zero out contributions in dead cells and core boundary
      do is = 0, ns-1
        call b2xzdd (nCv, 1, mpg, srw%sna0(1,0,is))
        call b2xzdd (nCv, 3, mpg, srw%smo0(1,0,is))
      enddo
      call b2xzdd (nCv, 1, mpg, srw%sne0)
      call b2xzdd (nCv, 3, mpg, srw%she0)
      call b2xzdd (nCv, 3, mpg, srw%shi0)
      call b2xzdd (nCv, 3, mpg, srw%shn0)
      call b2xzdd (nCv, 3, mpg, srw%sch0)
*   ..add to previous contributions
      call b2saxpy (nCv*2*ns, 1.0_R8, srw%sna0, 1, sr%sna, 1)
      call b2saxpy (nCv*4*ns, 1.0_R8, srw%smo0, 1, sr%smo, 1)
      call b2saxpy (nCv*2, 1.0_R8, srw%sne0, 1, sr%sne, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%she0, 1, sr%she, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%shi0, 1, sr%shi, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%shn0, 1, sr%shn, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%sch0, 1, sr%sch, 1)



*   ..compute sources by inputfiles (heating,...)
      if (sources_inputfile.eq.1) then
        call source_input(nx,ny,ns,
     &   sna0,smo0,she0,shi0,sch0,sne0,vol,ncall)
*   ..zero out contributions in dead cells and core boundary
        do is = 0, ns-1
          call b2xzdd (nCv, 1, mpg, srw%sna0(1,0,is))
          call b2xzdd (nCv, 3, mpg, srw%smo0(1,0,is))
        enddo
        call b2xzdd (nCv, 1, mpg, srw%sne0)
        call b2xzdd (nCv, 3, mpg, srw%she0)
        call b2xzdd (nCv, 3, mpg, srw%shi0)
        call b2xzdd (nCv, 3, mpg, srw%sch0)
*   ..add to previous contributions
        call b2saxpy (nCv*2*ns, 1.0_R8, srw%sna0, 1, sr%sna, 1)
        call b2saxpy (nCv*4*ns, 1.0_R8, srw%smo0, 1, sr%smo, 1)
        call b2saxpy (nCv*2, 1.0_R8, srw%sne0, 1, sr%sne, 1)
        call b2saxpy (nCv*4, 1.0_R8, srw%she0, 1, sr%she, 1)
        call b2saxpy (nCv*4, 1.0_R8, srw%shi0, 1, sr%shi, 1)
        call b2saxpy (nCv*4, 1.0_R8, srw%sch0, 1, sr%sch, 1)
      endif
#endif


*   ..compute sources due to electron-atom processes
*     (parallel momentum source must be pure friction)
      call b2stel (nCv, nFc, nVx, ns, ismain, BoRiS,
     &  geo, mpg, pl, dv, rt, rtw, srw) !srv 11.09.09 {
*   ..zero out contributions in dead cells and core boundary
      do is = 0, ns-1
        call b2xzdd (nCv, 1, mpg, srw%sna0(1,0,is))
        call b2xzdd (nCv, 3, mpg, srw%smq0(1,0,is))
      enddo
      call b2xzdd (nCv, 3, mpg, srw%she0)
      call b2xzdd (nCv, 3, mpg, srw%shi0)
      call b2xzdd (nCv, 3, mpg, srw%shn0)
*   ..add to previous contributions
      call b2saxpy (nCv*2*ns, 1.0_R8, srw%sna0, 1, sr%sna, 1)
      call b2saxpy (nCv*4*ns, 1.0_R8, srw%smq0, 1, sr%smq, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%she0, 1, sr%she, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%shi0, 1, sr%shi, 1)
      call b2saxpy (nCv*4, 1.0_R8, srw%shn0, 1, sr%shn, 1)
*     (no contribution to sch or sne)
*   ..compute sources due to charge exchange
      rcxnareg=0.0_R8
      rcxhireg=0.0_R8
      rcxmoreg=0.0_R8


!WG_TODO      rcxna=0.0_R8
!WG_TODO      rcxhi=0.0_R8
!WG_TODO      rcxmo=0.0_R8

      if (iout_b2wdat .eq. 4) then  ! IYS 06.11.2017
        if (.not. allocated(shi0_cx)) then
          write(*,*) 'Allocating arrays for cx'
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
        call b2sqcx (nCv, ns, ev, am(iscx(k)), pl%ti, pl%tn,
     &   rt%rlcx(1,0,0,k), rtw%rcx)
*    ..compute sources due to charge exchange
*      (parallel momentum source must be pure friction)
        call b2stcx (nCv, nFc, ns, iscx(k), ismain, BoRiS, geo, mpg,
     &    pl%na, pl%ua, pl%ti, pl%tn,  dv%ni, rtw%rcx, srw%sna0, 
     &    srw%smq0, srw%shi0, srw%shn0)       !srv 07.07.09
*   ..zero out contributions in dead cells and core boundary
        do is = 0, ns-1
         call b2xzdd (nCv, 1, mpg, srw%sna0(1,0,is))
         call b2xzdd (nCv, 3, mpg, srw%smq0(1,0,is))
        enddo
        call b2xzdd (nCv, 3, mpg, srw%shi0)
        call b2xzdd (nCv, 3, mpg, srw%shn0)
*    ..add to previous contributions
        call b2saxpy (nCv*2*ns, 1.0_R8, srw%sna0, 1, sr%sna, 1)
        call b2saxpy (nCv*4*ns, 1.0_R8, srw%smq0, 1, sr%smq, 1)
        call b2saxpy (nCv*4, 1.0_R8, srw%shi0, 1, sr%shi, 1)
        call b2saxpy (nCv*4, 1.0_R8, srw%shn0, 1, sr%shn, 1)
*      (no contribution to she, sne or sch)
        if (iout_b2wdat .eq. 4) then  ! IYS 06.11.2017
          write(*,*) 'Summing cx sources over hydrogenic neutrals'
          call b2saxpy (nCv*2*ns, 1.0_R8, srw%sna0, 1, sna0_cx, 1)
          call b2saxpy (nCv*4*ns, 1.0_R8, srw%smq0, 1, smo0_cx, 1)
          call b2saxpy (nCv*4, 1.0_R8, srw%shi0, 1, shi0_cx, 1)
        endif                         ! IYS 06.11.2017
      enddo


      if (iout_b2wdat .eq. 4) then  ! IYS 06.11.2017
        write(*,*) 'Writing cx sources'
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          wrk(:) = sna0_cx(:,0,is)+sna0_cx(:,1,is)*pl%na(:,is)
          call my_out_us(70,nCv,0,wrk,'b2stcx_sna_'//chns)
          wrk(:) = shi0_cx(:,0)+shi0_cx(:,1)*pl%ti+
     &             shi0_cx(:,2)*dv%ni(:,0)+
     &             shi0_cx(:,3)*pl%ti*dv%ni(:,0)
          call my_out_us(70,nCv,0,wrk,'b2stcx_shi')
          wrk(:) = smo0_cx(:,0,is)+
     &             smo0_cx(:,1,is)*pl%ua(:,is)+
     &             smo0_cx(:,2,is)*pl%na(:,is)*mp*am(is)+
     &             smo0_cx(:,3,is)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
	  call my_out_us(70,nCv,0,wrk,'b2stcx_smq'//chns)
        enddo
        if (allocated(shi0_cx)) then
          write(*,*) 'Deallocating arrays for cx'
          deallocate (shi0_cx)
          deallocate (sna0_cx)
          deallocate (smo0_cx)
        endif
      endif

      if(iout.ne.0) then                                                 !srv 30.06.08 {
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          do k=0,1
            write (chk,'(i1)') k
            call my_out_us(70,nCv,0,sr%sna(1,k,is),
     &       'b2sral_sna'//chk//chns)
          enddo
          do k=0,3
            write (chk,'(i1)') k
            call my_out_us(70,nCv,0,sr%smo(1,k,is),
     &       'b2sral_smo'//chk//chns)
            call my_out_us(70,nCv,0,sr%smq(1,k,is),
     &       'b2sral_smq'//chk//chns)
          enddo
        enddo
        do k=0,3
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,sr%she(1,k),'b2sral_she'//chk)
          call my_out_us(70,nCv,0,sr%shi(1,k),'b2sral_shi'//chk)
          call my_out_us(70,nCv,0,sr%shn(1,k),'b2sral_shn'//chk)
          call my_out_us(70,nCv,0,sr%sch(1,k),'b2sral_sch'//chk)
#ifdef WG_TODO
          if (use_astra.ne.0) then
            call my_out_us(70,nCv,0,sheNb_astra(1,k),
     &       'b2sral_sheNb_astra'//chk)
            call my_out_us(70,nCv,0,shiNb_astra(1,k),
     &       'b2sral_shiNb_astra'//chk)
          endif
#endif
        enddo
        do k=0,1
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,sr%sne(1,k),'b2sral_sne'//chk)
        enddo
      endif                                                              !srv 30.06.08 }

csw 17oct2012, erase sna, smo, smq, (she, shi, sch) for fluid neutrals
      if(no_fluid_neutral_sources .eq. 1) then
       do is=0,ns-1
        if(is_neutral(is)) then
          sr%sna(:,:,is) =0.d0
          sr%smo(:,:,is) =0.d0
          sr%smq(:,:,is) =0.d0
        endif
       enddo
      endif

* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sral

      end subroutine b2sral
