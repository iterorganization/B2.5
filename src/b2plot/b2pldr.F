      subroutine b2pldr (ninp, nout, ns)
      use b2mod_types
      use b2mod_version
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmgs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_b2cmrc
      use b2mod_b2cmwg
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
*   ..read geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
      call read_b2fgmtry (ninp(1), mpg, geo)
      call init_geometry (mpg, geo, switch)
*   ..read physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, ns)
*   ..read revised parameters in physics common
*     (ninp(0) need supply only parameters that are to be changed from
*     the values read in as default physics common.)
      call b2rfcp (ninp(0), ns)
*   ..test parameters
      call b2xvcp (ns)
*   ..read and echo code internal parameters
    1 continue
       read (ninp(0),'(a)',end=2,err=91) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         call ipsetc (cnamip, cvalip)
       end if
       goto 1
    2 continue
      call read_switch_afn(switch)
      call set_defaults_switches(switch,ns)
      call read_switches(switch)
      call check_switches(switch,ns)
*   ..read plasma state file
      call cfruch (ninp(3), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(3), b2fstate_version,
     &  ns, zamin, zamax, zn, am, .true.)
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
