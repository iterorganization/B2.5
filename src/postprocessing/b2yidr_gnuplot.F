*-----------------------------------------------------------------------
*.specification

      subroutine b2yidr_gnuplot (ninp, nout, nx, ny, ns)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_rates
      use b2mod_ppout
      use b2mod_version
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_work
      use b2mod_b2cmpa
      use b2mod_b2cmrc
      use b2mod_b2cmfs
      use b2mod_ppout
      use b2mod_geo_corner
      use b2mod_work
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
*     B2YIDR is a post-processor for b2. It reads and displays the
*     plasma state. This routine does not make use of any metric
*     information.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine is a post-processor for the b2 code. It reads the
*     result of a b2 calculation, and produces printed and graphical
*     summaries. This routine does not require any metric information.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:4) integer array, input.
*     ninp specifies the input unit numbers as follows.
*     ninp(0): formatted; provides output control parameters.
*     ninp(1): un*formatted; provides the plasma state.
*     ninp(2): un*formatted; provides atomic rate data.
*     ninp(3): un*formatted; provides the geometry.
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
      integer is, idum(0:9), ismain
      real (kind=R8) :: t0, t1, 
     *  rsa(-1:nx,-1:ny,0:ns-1), rra(-1:nx,-1:ny,0:ns-1),
     &  rqa(-1:nx,-1:ny,0:ns-1), rrd(-1:nx,-1:ny,0:ns-1),
     &  rbr(-1:nx,-1:ny,0:ns-1)
      character lblmn*120, lblrc*120, lblgm*120
      character cnamip*80, cvalip*80
*   ..initialisation
      save ismain
      data ismain/1/
*   ..procedures
      intrinsic abs
      logical ltst
      external subini, subend, xertst, cfruch, b2rugm, cfruin, b2xvps
      external b2rups, b2ruzd, b2wgpi_gnuplot, b2rurc, b2spel, b2sqel
      external ipsetc, ipgeti, b2rfcp, b2rflb, b2xxid
      ltst(t0,t1) = abs(t0-t1).le.1.0e-6_R8*(abs(t0)+abs(t1))

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
*     na - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), na(ix,iy,is) specifies
*     the density of atomic species (is) on the (ix,iy) cell.
*     It will hold that 0.lt.na(,,).
*
*     ne - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ne(ix,iy) specifies
*     the electron density on the (ix,iy) cell.
*     It will hold that 0.lt.ne(,).
*
*     ua - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), ua(ix,iy,is) specifies
*     the parallel velocity of atomic species (is) on the (ix,iy) cell.
*
*     uadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), uadia(ix,iy,0,is) specifies
*     the diamagnetic drift velocity of species (is) on the (ix,iy) cell, and
*     uadia(ix,iy,1,is) the radial drift velocity.
* 
*     te - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), te(ix,iy) specifies the electron
*     temperature on the (ix,iy) cell. It will hold that 0.lt.te(,).
*
*     ti - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), ti(ix,iy) specifies the all atom
*     temperature on the (ix,iy) cell. It will hold that 0.lt.ti(,).
*
*     po - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), po(ix,iy) specifies the electric
*     potential on the (ix,iy) cell.
*
*     fna - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let fnax(,,)=fna(,,0,) and fnay(,,)=fna(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), fnax(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fnax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fnay(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fnay(-1:nx,iy,0:ns-1) will hold 0.
*
*     fhe - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhex(,)=fhe(,,0) and fhey(,)=fhe(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhex(ix,iy) specifies the electron
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhex(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhey(ix,iy) specifies the electron
*     heat flux through the face between the (ix,iy-1) and (ix,iy)
*     cells. For iy.eq.-1, fhey(-1:nx,iy) will hold 0.
*
*     fhi - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhix(,)=fhi(,,0) and fhiy(,)=fhi(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhix(ix,iy) specifies the all atom
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhiy(ix,iy) specifies the all atom
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For iy.eq.-1, fhiy(-1:nx,iy) will hold 0.
*
*     fch - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchx(,)=fch(,,0) and fchy(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchx(ix,iy) specifies the electric
*     current through the face between the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchy(ix,iy) specifies the electric
*     current through the face between the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fchy(-1:nx,iy) will hold 0.
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
      call subini ('b2yidr_gnuplot')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.
     &  1.le.ninp(3).and.1.le.ninp(4).and.
     &  1.le.nout(0).and.1.le.nout(1), 
     &  'faulty argument ninp, nout')
*   ..test dimensions
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

* ..computation
*   ..obtain lblmn
      call b2xxid ('b2mn', lblmn(1:60))
      call b2rflb (ninp(4), lblmn(61:120))
*   ..read revised parameters in physics common
*     (ninp(4) need supply only parameters that are to be changed from
*     the values read in as default physics common.)
      call b2rfcp (ninp(4), ns)
*   ..read and echo code internal parameters
    1 continue
       read (ninp(4),*,end=2,err=91) cnamip, cvalip
       call ipsetc (cnamip, cvalip)
       goto 1
    2 continue
*   ..read plasma state file
      call cfruch (ninp(1), 120, lblmn, 'label')
*   ..read header for geometry
      call cfruch (ninp(3), 120, lblgm, 'label')
      call cfruin (ninp(3), 1, idum, 'isymm')
      isymm = idum(0)
*   ..read geometry and metric
      call b2rugm (ninp(3), nx, ny, crx, cry, fpsi,
     &  ffbz, bb, vol, hx, hy, qz, qc, gs, pbs, wbbl, wbbr, wbbv)
*   ..read species data
      call b2ruzd (ninp(1), b2fstate_version, ns, zamin, zamax, zn, am)
*   ..read plasma state
      call b2rups (ninp(1), nx, ny, ns, ne, na, ua, uadia, te,		 !xpb
     &  ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     1	b2fstate_version)
*   ..test plasma state and fluxes
      call b2xvps (nx, ny, ns, 
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
*   ..read atomic rate data
      call cfruch (ninp(2), 120, lblrc, 'label')
      call read_b2mod_ppout('rates', ninp(2), lblrc)
      call b2rurc (ninp(2), b2frates_version)
*    ..compute log-log electron rate coefficients
      call b2spel (nx, ny, ns, ev, te, ne, rlsa, rlra,
     &     rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi)
*   ..compute electron rate coefficients
      call ipgeti ('b2mndr_ismain', ismain)
      call xertst (0.le.ismain.and.ismain.lt.ns, 
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
      call b2sqel (nx, ny, ns, ismain, ev, te, 
     &    rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &    rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi, wrk0)
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

*   ..plot state
      call b2wgpi_gnuplot (nout(0), nx, ny, ns, ne, na, ua, te,
     &  ti, po, fna, fhe, fhi, fch, kinrgy,
     &  gs, qc, ni)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read internal parameters')

*-----------------------------------------------------------------------
*.end b2yidr_gnuplot

      end subroutine b2yidr_gnuplot
