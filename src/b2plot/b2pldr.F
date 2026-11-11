      subroutine b2pldr (ninp, nout, ns)
      use b2mod_types
      use b2mod_version
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_b2cmrc
      use b2mod_ppout
      use b2mod_diag
      use b2mod_math
      use b2mod_ipmain
      use b2mod_subsys
      use b2mod_switches
      use b2us_io
      use b2us_geo
      use b2us_map
      use b2us_data
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:6), nout(0:2), ns
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*   ..local variables
      integer is, idum(0:9), k, ismain
      integer external_species
      character lblgm*120, lblcp*120, lblmn*120, lblrc*120
      character cnamip*80, cvalip*80, line*500
*   ..initialisation
      save ismain, external_species
      data ismain/1/, external_species/0/
*   ..procedures
      logical is_comment, cfendf
      external xertst, xerrab, is_comment
      external cfruch, cfruin, b2rfcp, ipgeti, cfendf
      external b2rups, b2rugm, b2rucp, b2plot, b2xxid, b2rflb, b2rurc
      external b2ruzd, b2spel, b2sqel, b2sqcx, b2spcx, b2xvps, b2xvcp

*   ..subprogram start-up calls
      call subini ('b2pldr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.
     &             1.le.ninp(3).and.1.le.ninp(4).and.1.le.ninp(5).and.
     &             1.le.ninp(6).and.
     &             1.le.nout(0).and.1.le.nout(1).and.1.le.nout(2),
     &  'faulty argument ninp, nout')
*   ..test dimensions
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

* ..computation
*   ..obtain lblmn
      call b2xxid ('b2mn', lblmn(1:60))
      call b2rflb (ninp(0), lblmn(61:120))
      label = lblmn
*   ..read physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, ns)
*   ..read revised parameters in physics common
*     (ninp(0) need supply only parameters that are to be changed from
*     the values read in as default physics common.)
      call b2rfcp (ninp(0), ns)
*   ..test parameters
      call b2xvcp (ns)
*   ..read code internal parameters
    1 continue
       read (ninp(0),'(a)',end=2,err=91) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         call ipsetc (cnamip, cvalip)
       end if
       goto 1
    2 continue
*   ..assign switches
      call read_first_switches(switch)
      call set_defaults_switches(switch,ns)
      call read_switches(switch)
*   ..read geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
      call read_b2fgmtry (ninp(1), mpg, geo)
      call init_geometry (mpg, geo, switch)
*   ..read plasma state file
      call cfruch (ninp(3), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(3), b2fstate_version,
     &  ns, zamin, zamax, zn, am, .true.)
      call check_switches(switch,ns,nspecies)
*   ..read plasma state
      call read_b2fstate (ninp(3), mpg%nCv, mpg%nFc, ns, state)
*   ..test plasma state and fluxes
      call b2xvps (mpg%nCv, mpg%nFc, ns, state%pl, state%dv)
*   ..read atomic rate data
      call cfverr (ninp(5), b2frates_version)
      call cfruch (ninp(5), 120, lblrc, 'label')
      call read_b2mod_ppout('rates', ninp(5), lblrc)
      call b2rurc (ninp(5), b2frates_version)
*   ..test rtns, rtzmin, rtzmax, rtzn
      call xertst (ns.eq.rtns, 'faulty input rtns--wrong rate table?')
      do is = 0, ns-1
       call xertst (ltst(zamin(is),rtzmin(is)),
     &   'faulty input rtzmin--wrong rate table?')
       call xertst (ltst(zamax(is),rtzmax(is)),
     &   'faulty input rtzmax--wrong rate table?')
       call xertst (ltst(zn(is),rtzn(is)),
     &   'faulty input rtzn--wrong rate table?')
      enddo
*   ..read parameters from tracing file
      call cfruch (ninp(6), 120, lblmn, 'label')
      call b2ruzd (ninp(6), b2ftrack_version,
     &  ns, zamin, zamax, zn, am, .true.)
*   ..tracing loop
   10 if (.not.cfendf(ninp(6))) then
*    ..get itim, itcnt
       call cfruin (ninp(6), 2, idum, 'itim,itcnt')
*    ..get residuals and corrections
       call cfrure (ninp(6), mpg%nCv*ns, state%dv%resco, 'resco')
       call cfrure (ninp(6), mpg%nCv*ns, state%dv%corpa, 'corpa')
       call cfrure (ninp(6), mpg%nCv*ns, state%dv%resmo, 'resmo')
       call cfrure (ninp(6), mpg%nCv*ns, state%dv%corua, 'corua')
       call cfrure (ninp(6), mpg%nCv, state%dv%resmt, 'resmt')
       call cfrure (ninp(6), mpg%nCv, state%dv%corut, 'corut')
       call cfrure (ninp(6), mpg%nCv, state%dv%reshe, 'reshe')
       call cfrure (ninp(6), mpg%nCv, state%dv%corte, 'corte')
       call cfrure (ninp(6), mpg%nCv, state%dv%reshi, 'reshi')
       call cfrure (ninp(6), mpg%nCv, state%dv%corti, 'corti')
       call cfrure (ninp(6), mpg%nCv, state%dv%respo, 'respo')
       call cfrure (ninp(6), mpg%nCv, state%dv%corpo, 'corpo')
*    ..tracing loop
       goto 10
      endif
*    ..compute log-log electron rate coefficients
      call b2spel (mpg%nCv, ns, ev, state%pl%te, state%dv%ne, state%rt)
*   ..compute electron rate coefficients
      call ipgeti ('b2mndr_ismain', ismain)
      call xertst (0.le.ismain.and.ismain.lt.ns,
     &         'faulty b2mndr_ismain')
      call xertst (.not.is_neutral(ismain),
     &         'ismain cannot refer to a neutral species')
      call b2sqel (mpg%nCv, ns, ismain, switch, ev, state%pl%te,
     .      state%rt, state%rtw)
*   ..find nscx, iscx
*     (indices for neutral hydrogen species)
      state%rt%nscx = 0
      do is = 0, ns-1
        if (is_neutral(is).and.nint(zn(is)).eq.1) then
          call xertst (state%rt%nscx.lt.state%rt%nscxmax,
     1      'too many neutral hydrogenic species')
          state%rt%iscx(state%rt%nscx) = is
          state%rt%nscx = state%rt%nscx+1
        endif
      enddo
      do k = state%rt%nscx, state%rt%nscxmax-1
        state%rt%iscx(k) = -1
      enddo
*    ..compute log-log charge exchange rate coefficients
      do k = 0, state%rt%nscx-1
        call b2spcx (mpg%nCv, ns, ev, am(state%rt%iscx(k)),
     &   state%pl%ti, state%dv%ne,
     &   state%rt%rlcx(1:mpg%nCv,0:1,0:ns-1,k))
      enddo
*   ..compute charge exchange rate coefficients
      state%rtw%rcx = 0.0_R8
      do k = 0, state%rt%nscx-1
        call b2sqcx (mpg%nCv, ns, ev, switch, am(state%rt%iscx(k)),
     &   state%pl%ti, state%pl%tn,
     &   state%rt%rlcx(1:mpg%nCv,0:1,0:ns-1,k), state%rtw%rcx)
      enddo
*   ..obtain external species data if any
      call ipgeti ('b2mndr_external_species', external_species)
#ifdef WG_TODO
      if (external_species.ne.0) then
        call read_b2mod_external(nx,ny,ns,ismain)
      else
        call alloc_b2mod_external_dummy(nx,ny)
      endif
#endif
*   ..process
      call b2plot(ns)

* ..return
      call dealloc_b2mod_b2cmrc
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read internal parameters')

*-----------------------------------------------------------------------
*.end b2pldr

      end subroutine b2pldr

*-----------------------------------------------------------------------

      subroutine b2plrt (ns, switch, changed)
      use b2mod_types
      use b2mod_b2cmrc
      use b2mod_b2cmpa
      use b2mod_switches
      use b2mod_neutrals_namelist
      implicit none
*
*     The changes b2mn makes to the rate tables of b2frates when it
*     starts (b2mod_driver, b2mndr_0), so that b2plot computes the
*     rates the run used: the coronal model (b2mndr_coronal_model),
*     b2.atomic_physics_rescale.parameters (b2mndr_atomic_physics_
*     rescale), the recombination to the neutral species in runs
*     coupled to Eirene, which does that volume recombination
*     (b2mndr_rescale_neutrals_sources, times strasclfl of its 'V'
*     stratum), and hot electrons (b2mndr_dpc_mod_rates_te_hot,
*     b2mndr_dpc_mod_rates_ne_hot_frac). Keep it as there.
*     changed: a table was changed (the rates must be computed again).
*
      integer ns
      type (switches), intent(in) :: switch
      logical changed
      integer coronal, rescale_flag, is, istra
      real (kind=R8) :: scl, te_hot, ne_hot_frac
      real (kind=R8) ::
     1  rescale_sa(0:nsdecl-1), rescale_ra(0:nsdecl-1),
     2  rescale_qa(0:nsdecl-1), rescale_cx(0:nsdecl-1),
     3  rescale_rd(0:nsdecl-1), rescale_br(0:nsdecl-1)
      character*256 filename
      logical ex, streql
      external ipgeti, ipgetr, xertst, find_file, streql
      namelist /atomic_physics_rescale/
     & rescale_sa, rescale_ra, rescale_qa, rescale_cx,
     & rescale_rd, rescale_br

      changed = .false.
*   ..coronal model
      coronal = 0
      call ipgeti ('b2mndr_coronal_model', coronal)
      if (coronal.ne.0) then
        write(*,*)
     1   'Modifying rate coefficients to suppress density dependence'
        do is = 1, rtnn
          rtlsa(:,is,:)=rtlsa(:,0,:)
          rtlra(:,is,:)=rtlra(:,0,:)
          rtlqa(:,is,:)=rtlqa(:,0,:)
          rtlqr(:,is,:)=rtlqr(:,0,:)
          rtlcx(:,is,:)=rtlcx(:,0,:)
          rtlrd(:,is,:)=rtlrd(:,0,:)
          rtlbr(:,is,:)=rtlbr(:,0,:)
          rtlza(:,is,:)=rtlza(:,0,:)
          rtlz2(:,is,:)=rtlz2(:,0,:)
          rtlpt(:,is,:)=rtlpt(:,0,:)
          rtlpi(:,is,:)=rtlpi(:,0,:)
        enddo
        changed = .true.
      endif
*   ..rescaled atomic physics
      rescale_flag = 0
      call ipgeti ('b2mndr_atomic_physics_rescale', rescale_flag)
      if (rescale_flag.ne.0) then
        rescale_sa = 1.0_R8
        rescale_ra = 1.0_R8
        rescale_qa = 1.0_R8
        rescale_cx = 1.0_R8
        rescale_rd = 1.0_R8
        rescale_br = 1.0_R8
        filename = 'b2.atomic_physics_rescale.parameters'
        call find_file(filename,ex)
        open(99,file=filename)
        write(*,*) 'Reading ',trim(filename)
        read(99,atomic_physics_rescale)
        close(99)
        call xertst (min(
     1       minval(rescale_sa(0:rtns-1)),
     2       minval(rescale_ra(0:rtns-1)),
     3       minval(rescale_qa(0:rtns-1)),
     4       minval(rescale_cx(0:rtns-1)),
     5       minval(rescale_rd(0:rtns-1)),
     6       minval(rescale_br(0:rtns-1)))
     7       .gt.0.0_R8,
     8  'faulty value set in b2.atomic_physics_rescale.parameters')
        do is = 0, rtns-1
          rtlsa(:,:,is)=rtlsa(:,:,is)+log(rescale_sa(is))
          rtlra(:,:,is)=rtlra(:,:,is)+log(rescale_ra(is))
          rtlqa(:,:,is)=rtlqa(:,:,is)+log(rescale_qa(is))
          rtlqr(:,:,is)=rtlqr(:,:,is)+log(rescale_qa(is))
          rtlcx(:,:,is)=rtlcx(:,:,is)+log(rescale_cx(is))
          rtlrd(:,:,is)=rtlrd(:,:,is)+log(rescale_rd(is))
          rtlbr(:,:,is)=rtlbr(:,:,is)+log(rescale_br(is))
        enddo
        changed = .true.
      endif
*   ..recombination to the neutral species: done by Eirene
      if (switch%use_eirene.ne.0 .and.
     &    switch%neutral_sources_rescale.ne.1.0_R8) then
        write(*,*) 'Adding ',log(switch%neutral_sources_rescale),
     &   ' to relevant rtlra and rtlqr'
        do is = 0, ns-2
          scl = switch%neutral_sources_rescale
          do istra = 1, nstrai
            if (streql(crcstra(istra),'V').and.
     &        (species_start(istra).eq.is)) then
              scl = scl*strasclfl(istra)
            endif
          enddo
          if (is_neutral(is) .and. lnext(is,is+1)) then
            rtlra(:,:,is+1)=rtlra(:,:,is+1)+log(scl)
            rtlqr(:,:,is+1)=rtlqr(:,:,is+1)+log(scl)
          endif
        enddo
        changed = .true.
      endif
*   ..hot electrons
      te_hot = 0.0_R8
      ne_hot_frac = 0.0_R8
      call ipgetr ('b2mndr_dpc_mod_rates_te_hot', te_hot)
      call ipgetr ('b2mndr_dpc_mod_rates_ne_hot_frac', ne_hot_frac)
      if (te_hot.gt.0.0_R8.and.ne_hot_frac.gt.0.0_R8) then
        call b2plrt_hot (te_hot, ne_hot_frac)
        changed = .true.
      endif
      return
      end subroutine b2plrt

*-----------------------------------------------------------------------
*     As dpc_mod_rates of b2mod_driver: rates with a fraction
*     ne_hot_frac of electrons at te_hot.

      subroutine b2plrt_hot(te_hot, ne_hot_frac)
      use b2mod_types
      use b2mod_b2cmrc

      implicit none

      real (R8) :: te_hot
      real (R8) :: ne_hot_frac
      real (R8) :: rlte, fxte
      real (R8) :: rate_sa, rate_ra, rate_qa, rate_rd, rate_br,
     *             rate_za, rate_z2, rate_pt, rate_pi, rate_qr
!!    real (R8) :: rate_cx
      integer ite, ine, ise, ixte

      write(*,*) 'b2plrt_hot: ', te_hot, ne_hot_frac

      ite = 0
      rlte = max(rtlt(0),min(rtlt(rtnt),log(te_hot)))

    1 if (ite.lt.rtnt-1) then
        if (rtlt(ite+1).le.rlte) then
          ite = ite+1
          goto 1
        endif
      endif
    2 if (0.lt.ite) then
        if (rlte.lt.rtlt(ite)) then
          ite = ite-1
          goto 2
        endif
      endif
*     ..compute coefficient for linear interpolation
      ixte = ite
      fxte = (rlte-rtlt(ite))/(rtlt(ite+1)-rtlt(ite))
      write(*,*) ixte
      write(*,*) rtn(rtnn/2)
      write(*,*) rtt(:)
      write(*,*) exp(rtlsa(:,rtnn/2,0))
      do ise= 0, rtns-1
        do ine= 0, rtnn
          rate_sa = exp((1.0_R8-fxte)*rtlsa(ixte,ine,ise) +
     1     fxte*rtlsa(ixte+1,ine,ise)) * ne_hot_frac
          rate_ra = exp((1.0_R8-fxte)*rtlra(ixte,ine,ise) +
     1     fxte*rtlra(ixte+1,ine,ise)) * ne_hot_frac
          rate_qa = exp((1.0_R8-fxte)*rtlqa(ixte,ine,ise) +
     1     fxte*rtlqa(ixte+1,ine,ise)) * ne_hot_frac
          rate_qr = exp((1.0_R8-fxte)*rtlqr(ixte,ine,ise) +
     1     fxte*rtlqr(ixte+1,ine,ise)) * ne_hot_frac
          rate_rd = exp((1.0_R8-fxte)*rtlrd(ixte,ine,ise) +
     1     fxte*rtlrd(ixte+1,ine,ise)) * ne_hot_frac
          rate_br = exp((1.0_R8-fxte)*rtlbr(ixte,ine,ise) +
     1     fxte*rtlbr(ixte+1,ine,ise)) * ne_hot_frac
!!        rate_cx = exp((1.0_R8-fxte)*rtlcx(ixte,ine,ise) +
!!   1     fxte*rtlcx(ixte+1,ine,ise)) * ne_hot_frac
          rate_za = ((1.0_R8-fxte)*rtlza(ixte,ine,ise) +
     1     fxte*rtlza(ixte+1,ine,ise)) * ne_hot_frac
          rate_z2 = ((1.0_R8-fxte)*rtlz2(ixte,ine,ise) +
     1     fxte*rtlz2(ixte+1,ine,ise)) * ne_hot_frac
          rate_pt = ((1.0_R8-fxte)*rtlpt(ixte,ine,ise) +
     1     fxte*rtlpt(ixte+1,ine,ise)) * ne_hot_frac
          rate_pi = ((1.0_R8-fxte)*rtlpi(ixte,ine,ise) +
     1     fxte*rtlpi(ixte+1,ine,ise)) * ne_hot_frac
          do ite= 0, rtnt
            rtlsa(ite,ine,ise) =
     1       log(exp(rtlsa(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_sa)
            rtlra(ite,ine,ise) =
     1       log(exp(rtlra(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_ra)
            rtlqa(ite,ine,ise) =
     1       log(exp(rtlqa(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_qa)
            rtlqr(ite,ine,ise) =
     1       log(exp(rtlqr(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_qr)
            rtlrd(ite,ine,ise) =
     1       log(exp(rtlrd(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_rd)
            rtlbr(ite,ine,ise) =
     1       log(exp(rtlbr(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_br)
!!          rtlcx(ite,ine,ise) =
!!   1       log(exp(rtlcx(ite,ine,ise))*(1.0_R8-ne_hot_frac)+rate_cx)
            rtlza(ite,ine,ise) =
     1       rtlza(ite,ine,ise)*(1.0_R8-ne_hot_frac)+rate_za
            rtlz2(ite,ine,ise) =
     1       rtlz2(ite,ine,ise)*(1.0_R8-ne_hot_frac)+rate_z2
            rtlpt(ite,ine,ise) =
     1       rtlpt(ite,ine,ise)*(1.0_R8-ne_hot_frac)+rate_pt
            rtlpi(ite,ine,ise) =
     1       rtlpi(ite,ine,ise)*(1.0_R8-ne_hot_frac)+rate_pi
          enddo
        enddo
      enddo
      write(*,*) exp(rtlsa(:,rtnn/2,0))
      return
      end subroutine b2plrt_hot
