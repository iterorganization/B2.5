*-----------------------------------------------------------------------
*.specification

      subroutine b2yndr (ninp, nout, nx, ny, ns)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_work
      use b2mod_indirect
      use b2mod_version
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmrc
      use b2mod_ppout
      use b2mod_ipmain
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:4), nout(0:1), nx, ny, ns
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YNDR is a postprocessor for B2. It displays the progress of
*     the inner iterations.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine is a postprocessor for the B2 code. It displays the
*     evolution of the corrections and residuals. This routine does not
*     require any metric information.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:4) integer array, input.
*     ninp specifies the input unit numbers as follows.
*     ninp(0): formatted; provides output control parameters.
*     ninp(1): un*formatted; provides the corrections and residuals.
*     ninp(2): un*formatted; provides atomic rate data.
*     ninp(3): un*formatted; provides the plasma state.
*     ninp(4): formatted; provides run switches.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides printed output.
*     nout(1): graphical; provides graphical output.
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
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer n2, is, itim, itcnt, idum(0:9), iftp, ismain
      real (kind=R8) ::
     *  rsa(-1:nx,-1:ny,0:ns-1), rra(-1:nx,-1:ny,0:ns-1),
     &  rqa(-1:nx,-1:ny,0:ns-1), rrd(-1:nx,-1:ny,0:ns-1),
     &  rbr(-1:nx,-1:ny,0:ns-1)
      character lblmn*120, lblrc*120, idout*80, id*8
      character cnamip*80, cvalip*80, line*500
*   ..initialisation
      save ismain
      data ismain/1/
*   ..procedures
      logical streql, cfendf, is_comment
      external xertst, xerrab, streql, gacwk,
     &  gdawk, cfruch, cfruin, cfrure, cfendf, lbfilx, ipgeti,
     &  is_comment
      external b2wgcr, b2ruzd, b2rurc, b2spel, b2sqel, b2rups, b2xvps,
     &  b2rfcp, b2rflb, b2xxid

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
*      me, mp, ev, qe, eps0, mu0 - real*8 constant.
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
*     resco - (-1:nx,-1:ny,0:ns-1) real*8 array, input.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), resco(ix,iy,is) specifies
*     the residual of the continuity equation on the (ix,iy) cell for
*     species (is).
*
*     resmo - (-1:nx,-1:ny,0:ns-1) real*8 array, input.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), resmo(ix,iy,is) specifies
*     the residual of the parallel momentum equation on the (ix,iy) cell
*     for species (is).
*
*     resmt - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), resmt(ix,iy) specifies the residual
*     of the total parallel momentum equation on the (ix,iy) cell.
*
*     reshe - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), reshe(ix,iy) specifies the residual
*     of the electron heat balance equation on the (ix,iy) cell.
*
*     reshi - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), reshi(ix,iy) specifies the residual
*     of the all atom heat balance equation on the (ix,iy) cell.
*
*     respo - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), respo(ix,iy) specifies the
*     residual of the electric potential equation on the (ix,iy) cell.
*
*     corpa - (-1:nx,-1:ny,0:ns-1) real*8 array, input.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), corpa(ix,iy,is) specifies
*     the last computed correction to the pressure of species (is) on
*     the (ix,iy) cell. corpa is computed by solving an equation that
*     has resco as its right hand side. It is not necessarily the same
*     as the pressure correction over one complete iteration.
*
*     corua - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), corua(ix,iy,is) specifies
*     the last computed correction to the parallel velocity ua of
*     species (is) on the (ix,iy) cell.
*     corua is computed by solving an equation that has resmo as its
*     right hand side. It is not necessarily the same as the change in
*     the parallel velocity over one complete iteration.
*
*     corut - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (0:nx,-1:ny), corut(ix,iy) specifies the last
*     computed correction to the parallel velocity ua of all species
*     simultaneously on the (ix,iy) cell.
*     corut is computed by solving an equation that has resmt as its
*     right hand side. It is not necessarily the same as the change in
*     the parallel velocity over one complete iteration.
*
*     corte - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), corte(ix,iy) specifies the last
*     computed correction to the electron temperature on the (ix,iy)
*     cell. corte is computed by solving an equation that has reshe as
*     its right hand side. It is not necessarily the same as the change
*     in the electron temperature over one complete iteration.
*
*     corti - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), corti(ix,iy) specifies the last
*     computed correction to the atom temperature on the (ix,iy) cell.
*     corti is computed by solving an equation that has reshi as its
*     right hand side. It is not necessarily the same as the change in
*     the atom temperature over one complete iteration.
*
*     corpo - (-1:nx,-1:ny) real*8 array, input.
*     For (ix,iy) in (-1:nx,-1:ny), corpo(ix,iy) specifies the last
*     computed correction to the electric potential on the (ix,iy)
*     cell. corpo is computed by solving an equation that has respo as
*     its right hand side. It is not necessarily the same as the change
*     in the electric potential over one complete iteration.
*
*     lblmn - character*120
*     lblmn identifies the instance of the program employed for the main
*     calculation.
*
*     ismain - integer.
*     ismain identifies the index of the main plasma species.
*     This cannot be a neutral species.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2yndr')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.
     &             1.le.ninp(3).and.1.le.ninp(4).and.
     &             1.le.nout(0).and.1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..test dimensions
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

* ..computation
*   ..read output control
      read (ninp(0),'(a8)',err=91) id
      call xertst (streql(id,'*idout'),
     &  'faulty input instead of *idout, id='//id)
      read (ninp(0),*,err=92) iftp, idout
*   ..obtain lblmn
      call b2xxid ('b2mn', lblmn(1:60))
      call b2rflb (ninp(4), lblmn(61:120))
*   ..read revised parameters in physics common
*     (ninp(4) need supply only parameters that are to be changed from
*     the values read in as default physics common.)
      call b2rfcp (ninp(4), ns)
*   ..read and echo code internal parameters
    1 continue
       read (ninp(4),'(a)',end=2,err=93) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         call ipsetc (cnamip, cvalip)
       end if
       goto 1
    2 continue
*   ..read parameters
      call cfruch (ninp(1), 120, lblmn, 'label')
      call b2ruzd (ninp(1), b2ftrack_version,
     &  ns, zamin, zamax, zn, am, .true.)
*   ..write label
      write (nout(0),'(/2x,a)') 'label for this calculation:'
      write (nout(0),'(/4x,a/4x,a)') lblmn(1:60), lblmn(61:120)
*   ..write output control
      write (nout(0),'(/2x,a,2x,i3/2x,a,2x,a)')
     &  'style of contour plots:', iftp,
     &  'idout:', idout
*   ..read atomic rate data
      call cfruch (ninp(2), 120, lblrc, 'label')
      call get_AM_flag (lblrc, b2frates_flag)
      call read_b2mod_ppout('rates', ninp(2), b2frates_flag)
      call b2rurc (ninp(2), b2frates_version)
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
*   ..read plasma state file label
      call cfruch (ninp(3), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(3), b2fstate_version,
     &  ns, zamin, zamax, zn, am, .true.)
*   ..read plasma state
      call b2rups (ninp(3), nx, ny, ns, ne, na, ua, uadia, te,
     &  ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc,
     1  b2fstate_version)
*   ..test plasma state and fluxes
      call b2xvps (nx, ny, ns,
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
*    ..compute log-log electron rate coefficients
      call b2spel (nx, ny, ns, ev, te, ne, rlsa, rlra,
     &  rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi)
*   ..compute electron rate coefficients
      call ipgeti ('b2mndr_ismain', ismain)
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
      call b2sqel (nx, ny, ns, ismain, ev, te,
     &  rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &  rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi)

* ..produce plots
*   ..activate graphics on nout(1)
      call gacwk (nout(1))
*   ..mark contour style
      call lbfilx (iftp)
*   ..main loop
   10 if (.not.cfendf(ninp(1))) then
*    ..get itim, itcnt
       call cfruin (ninp(1), 2, idum, 'itim,itcnt')
       itim = idum(0)
       itcnt = idum(1)
*    ..echo itim, itcnt
       write (nout(0),'(2x,a,i6,4x,a,i6)')
     &   'itim:', itim, 'itcnt:', itcnt
*    ..get residuals and corrections
       call cfrure (ninp(1), n2*ns, resco, 'resco')
       call cfrure (ninp(1), n2*ns, corpa, 'corpa')
       call cfrure (ninp(1), n2*ns, resmo, 'resmo')
       call cfrure (ninp(1), n2*ns, corua, 'corua')
       call cfrure (ninp(1), n2, resmt, 'resmt')
       call cfrure (ninp(1), n2, corut, 'corut')
       call cfrure (ninp(1), n2, reshe, 'reshe')
       call cfrure (ninp(1), n2, corte, 'corte')
       call cfrure (ninp(1), n2, reshi, 'reshi')
       call cfrure (ninp(1), n2, corti, 'corti')
       call cfrure (ninp(1), n2, respo, 'respo')
       call cfrure (ninp(1), n2, corpo, 'corpo')
*    ..plot corrections and residuals
       call b2wgcr (nout(0), idout, nx, ny, ns, resco, resmo,
     &   resmt, reshe, reshi, respo, corpa,
     &   corua, corut, corte, corti, corpo)
*    ..loop
       goto 10
      endif
*   ..deactivate graphics on nout(1)
      call gdawk (nout(1))

* ..return
      call dealloc_b2mod_b2cmrc
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read heading *idout')
   92 call xerrab ('error trying to read output control')
   93 call xerrab ('error trying to read internal parameters')

*-----------------------------------------------------------------------
*.end b2yndr

      end subroutine b2yndr
