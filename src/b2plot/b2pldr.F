      subroutine b2pldr (ninp, nout, nx, ny, ns)
      use b2mod_types
      use b2mod_version
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmfs
      use b2mod_b2cmgs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_b2cmrc
      use b2mod_b2cmwg
      use b2mod_ppout
      use b2mod_diag
      use b2mod_work
      use b2mod_external
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:6), nout(0:2), nx, ny, ns
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*   ..local variables
      integer is, idum(0:9), k, nscx, iscx(0:nscxmax-1), ismain
      integer external_species
      integer itim, itcnt, n2
      character lblgm*120, lblcp*120, lblmn*120, lblrc*120
      character cnamip*80, cvalip*80, line*500
      real (kind=R8) :: t0, t1
      real (kind=R8) ::
     *  rsa(-1:nx,-1:ny,0:ns-1), rra(-1:nx,-1:ny,0:ns-1),
     &  rqa(-1:nx,-1:ny,0:ns-1), rrd(-1:nx,-1:ny,0:ns-1),
     &  rbr(-1:nx,-1:ny,0:ns-1), rcx(-1:nx,-1:ny,0:ns-1)
*   ..initialisation
      save ismain, external_species
      data ismain/1/, external_species/0/
*   ..procedures
      logical ltst, is_comment, cfendf
      intrinsic sqrt
      external subini, subend, xertst, xerrab, is_comment
      external cfruch, cfruin, cfverr, b2rfcp, ipgeti, cfendf
      external b2rups, b2rugm, b2rucp, b2plot, b2xxid, b2rflb, b2rurc
      external b2ruzd, b2spel, b2sqel, b2sqcx, b2spcx, b2xvps, b2xvcp
      ltst(t0,t1) = abs(t0-t1).le.1.0e-6_R8*(abs(t0)+abs(t1))

*   ..subprogram start-up calls
      call subini ('b2pldr')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.
     &             1.le.ninp(3).and.1.le.ninp(4).and.1.le.ninp(5).and.
     &             1.le.ninp(6).and.
     &             1.le.nout(0).and.1.le.nout(1).and.1.le.nout(2),
     &  'faulty argument ninp, nout')
*   ..test dimensions
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
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
      call b2rugm (ninp(1), nx, ny, crx, cry, fpsi, ffbz,
     &  bb, vol, hx, hy, qz, qc, qs, gs, pbs, wbbl, wbbr, wbbv,
     &  wbbc, cell_width, cell_height, gmap)
*   ..read physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, nx, ns)
*   ..read revised parameters in physics common
*     (ninp(0) need supply only parameters that are to be changed from
*     the values read in as default physics common.)
      call b2rfcp (ninp(0), ns)
*   ..test parameters
      call b2xvcp (nx, ns)
*   ..read and echo code internal parameters
    1 continue
       read (ninp(0),'(a)',end=2,err=91) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         call ipsetc (cnamip, cvalip)
       end if
       goto 1
    2 continue
*   ..read plasma state file
      call cfruch (ninp(3), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(3), b2fstate_version,
     &  ns, zamin, zamax, zn, am, .true.)
*   ..read plasma state
      call b2rups (ninp(3), nx, ny, ns, ne, na, ua, uadia, te,         !xpb
     &  ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc,
     &  b2fstate_version)
*   ..test plasma state and fluxes
      call b2xvps (nx, ny, ns,
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
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
     & ns, zamin, zamax, zn, am, .true.)
*   ..tracing loop
   10 if (.not.cfendf(ninp(6))) then
*    ..get itim, itcnt
       call cfruin (ninp(6), 2, idum, 'itim,itcnt')
       itim = idum(0)
       itcnt = idum(1)
*    ..get residuals and corrections
       call cfrure (ninp(6), n2*ns, resco, 'resco')
       call cfrure (ninp(6), n2*ns, corpa, 'corpa')
       call cfrure (ninp(6), n2*ns, resmo, 'resmo')
       call cfrure (ninp(6), n2*ns, corua, 'corua')
       call cfrure (ninp(6), n2, resmt, 'resmt')
       call cfrure (ninp(6), n2, corut, 'corut')
       call cfrure (ninp(6), n2, reshe, 'reshe')
       call cfrure (ninp(6), n2, corte, 'corte')
       call cfrure (ninp(6), n2, reshi, 'reshi')
       call cfrure (ninp(6), n2, corti, 'corti')
       call cfrure (ninp(6), n2, respo, 'respo')
       call cfrure (ninp(6), n2, corpo, 'corpo')
*    ..tracing loop
       goto 10
      endif
*    ..compute log-log electron rate coefficients
      call b2spel (nx, ny, ns, ev, te, ne, rlsa, rlra,
     &     rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi)
*   ..compute electron rate coefficients
      call ipgeti ('b2mndr_ismain', ismain)
      call xertst (0.le.ismain.and.ismain.lt.ns,
     &         'faulty b2mndr_ismain')
      call xertst (.not.is_neutral(ismain),
     &         'ismain cannot refer to a neutral species')
      call b2sqel (nx, ny, ns, ismain, ev, te,
     &     rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &     rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi)
*   ..find nscx, iscx
*     (indices for neutral hydrogen species)
      nscx = 0
      do is = 0, ns-1
        if (is_neutral(is).and.zn(is).eq.1.0_R8) then
          call xertst (nscx.lt.nscxmax,
     1     'too many neutral hydrogen species')
          iscx(nscx) = is
          nscx = nscx+1
        endif
      enddo
      do k = nscx, nscxmax-1
        iscx(k) = -1
      enddo
*    ..compute log-log charge exchange rate coefficients
      do k = 0, nscx-1
        call b2spcx (nx, ny, ns, ev, am(iscx(k)), ti, ne,
     &   rlcx(-1,-1,0,0,k))
      enddo
*   ..compute charge exchange rate coefficients
      do k = 0, nscx-1
        call b2sqcx (nx, ny, ns, ev, am(iscx(k)), ti,
     &   rlcx(-1,-1,0,0,k), rcx)
      enddo
*  .. interpolate magnetic field to cell corners                         !xpb
c!!!This interpolation scheme does not guarantee that bb(3)=sqrt(bb(0:2)**2)!!!
      call my_intcor1(nx,ny,2,vol,periodic_bc,nncut,
     * leftix,leftiy,rightix,rightiy,
     * bottomix,bottomiy,topix,topiy,region,bb,wbbc)                     !xpb
*   ..compute extrapolated wbbc(:,:,3) component
      wbbc(:,:,3)=sqrt(wbbc(:,:,0)**2+wbbc(:,:,1)**2+wbbc(:,:,2)**2)
*   ..obtain external species data if any
      call ipgeti ('b2mndr_external_species', external_species)
      if (external_species.ne.0) then
        call read_b2mod_external(nx,ny,ns,ismain)
      else
        call alloc_b2mod_external_dummy(nx,ny)
      endif
*   ..process
      call b2plot(nx,ny,ns)

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
