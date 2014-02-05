*-----------------------------------------------------------------------
*.specification

      subroutine b2ahdr (ninp, nout, nx, ns)
      use b2mod_types
      use b2mod_version
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:0), nout(0:1), nx, ns
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AHDR prepares and writes to an external file a table of default
*     physics parameters for b2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine prepares default physics parameters for b2. The
*     parameters are written to an un*formatted file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:0) integer array, input.
*     ninp specifies the input unit numbers, as follows.
*     ninp(0): formatted; provides physics parameters.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides summary printed output.
*     nout(1): un*formatted; provides default physics parameters.
*
*     nx - integer, input.
*     nx specifies the number of poloidal cells in the calculation.
*     The cells are indexed by (-1:nx).
*     It will hold that 1.le.nx.
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
      integer idum(0:9)
      character lblcp*120
*   ..procedures
      external subini, subend, xertst, cfwuin, cfwuch, cfverw
      external b2xxid, b2rflb, b2rfcp, b2wfcp, b2wucp, b2xvcp, b2xvps

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
*     /b2cmpb/ specifies physics parameters for boundary conditions.
*     For this purpose the boundary is divided up into regions, and for
*     each region a set of coefficients defines the local conditions.
*     For detailed description of /b2cmpb/ see the code, especially the
*     routine b2stbc.
*      cbnrmx - integer constant.
*     cbnrmx is an upper bound for the number of boundary regions nreg,
*     (nreg.eq.cbnrso+cbnrno+cbnrwe+cbnrea) and is used in dimension
*     specifications.
*      cbirso, cbirno, cbirwe, cbirea - integer.
*     cbirso(-no,-we,-ea) is the index of the first region on the
*     south (north, west, east) boundary.
*      cbnrso, cbnrno, cbnrwe, cbnrea - integer.
*     cbnrso(-no,-we,-ea) is the number of regions on the south
*     (north, west, east) boundary.
*     It will hold that 1.le.cbnrso, 1.le.cbnrno, 1.le.cbnrwe, and
*     1.le.cbnrea; and 0.eq.cbirso, cbirso+cbnrso.eq.cbirno,
*     cbirno+cbnrno.eq.cbirwe, cbirwe+cbnrwe.eq.cbirea.
*      cbrbrk - (0:cbnrmx-1) real*8 array.
*     cbrbrk defines the boundaries between the boundary regions.
*      cbsna - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns and 0.le.ireg.lt.nreg, the eight coefficients
*     cbsna(0:7,is,ireg) specify standard form boundary conditions for
*     the continuity equation for species (is) on the (ireg) boundary
*     region.
*      cbsmo - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns and 0.le.ireg.lt.nreg, the eight coefficients
*     cbsmo(0:7,is,ireg) specify standard form boundary conditions for
*     the parallel momentum balance equation for species (is) on the
*     (ireg) boundary region.
*      cbshi - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns and 0.le.ireg.lt.nreg, the eight coefficients
*     cbshi(0:7,is,ireg) specify standard form boundary conditions for
*     the contribution to that atom heat equation by species (is) on
*     the (ireg) boundary region.
*      cbshe - (0:7,0:cbnrmx-1) real*8 array.
*     For 0.le.ireg.lt.nreg, the eight coefficients cbshe(0:7,ireg)
*     specify standard form boundary conditions for the electron heat
*     equation on the (ireg) boundary region.
*      cbsch - (0:7,0:cbnrmx-1) real*8 array.
*     For 0.le.ireg.lt.nreg, the eight coefficients cbsch(0:7,ireg)
*     specify standard form boundary conditions for the electric
*     potential equation on the (ireg) boundary region.
*      cbrec - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns-1 and 0.le.ireg.lt.nreg, the eight coefficients
*     cbrec(0:7,is,ireg) specify the recycling process for species (is)
*     on the boundary region ireg. The exact interpretation of these
*     coefficients is defined by the routine that computes the source
*     terms due to the recycling process (b2stbr). Routine b2stbr is
*     problem dependent.
*      cbmsa - (0:7,0:cbnrmx-1) real*8 array.
*      cbmsc - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.ireg.lt.nreg, the eight coefficients cbmsa(0:7,ireg)
*     specify miscellaneous processes on the boundary region ireg.
*     For 0.le.is.lt.ns-1 and 0.le.ireg.lt.nreg, the eight coefficients
*     cbmsc(0:7,is,ireg) specify miscellaneous processes for species
*     (is) on the boundary region ireg. The interpretation of these
*     coefficients is defined by the routine that computes the source
*     terms due to miscellaneous processes (b2stbm). Routine b2stbm is
*     application dependent.
*
*   ..transport coefficients, parameters
*
*     /b2cmpt/ - common block.
*     (Note: the integer constant nsdecl, which comes with /b2cmpa/,
*     is used instead of (ns) in dimension specifications.)
*      cfdf0 - (0:7,0:ns-1) real*8 array.
*     cfdf0 contains some auxiliary coefficients.  For 0.le.is.lt.ns,
*     the local expression
*       df0(,,is) = cfdf0(0,is)+cfdf0(1,is)/(ni(,,0)/1.0e20)+
*         cfdf0(2,is)*sqrt(ti(,)/(am(is)*mp))/
*          (cfdf0(3,is)*ni(,,0)+cfdf0(4,is)*ne(,))
*     is used in the definition of the anomalous transport coefficients.
*     df0 is meant to represent a diffusivity.  If cfdf0(2,is).eq.0 then
*     the corresponding term in df0 is understood as 0 independent of
*     cfdf0(3:4,is).
*      cfdna - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfdna(,is) are the coefficients in the
*     local expression
*       difna(,,is) = cfdna(0,is)+cfdna(1,is)/(ne(,)/1.0e20)+
*         cfdna(2,is)*df0(,)+
*         cfdna(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous diffusivity with respect to
*     grad.na; the transport equation for species (is) contains a term
*     difna(,,is)*'grad.na'.
*     If cfdna(7,is) is not zero, the dna0 array is flux-scaled with
*     a scaling factor of cfdna(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfdna(4:6,) are not used in the default version of the code.
*      cfdpa - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfdpa(,is) are the coefficients in the
*     local expression
*       difpa(,,is) = cfdpa(0,is)+cfdpa(1,is)/(ne(,)/1.0e20)+
*         cfdpa(2,is)*df0(,,is)+
*         cfdpa(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous diffusivity with respect to
*     grad.pa; the transport equation for species (is) contains a term
*     (1/(rza(,,is)*te(,)+ti(,)))*difpa(,,is)*'grad.pa'.
*     If cfdpa(7,is) is not zero, the dpa0 array is flux-scaled with
*     a scaling factor of cfdpa(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfdpa(4:6,) are not used in the default version of the code.
*      cfvla - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfvla(,is) are the coefficients in the
*     local expression
*       vlay(,,is) = cfvla(0,is)+cfvla(1,is)/(ne(,)/1.0e20)+
*         cfvla(2,is)*df0(,,is)+
*         cfvla(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous y-velocity of species (is); the
*     transport equation for species (is) contains a term
*     vlay(,,is)*na(,,is).
*     If cfvla(7,is) is not zero, the vla0 array is flux-scaled with
*     a scaling factor of cfvla(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfvla(4:6,) are not used in the default version of the code.
*      cfvsa - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfvsa(,is) are the coefficients in the
*     local expression
*       dfvsa(,,is) = cfvsa(0,is)+cfvsa(1,is)/(ne(,)/1.0e20)+
*         cfvsa(2,is)*df0(,,is)+
*         cfvsa(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     transverse viscosity, so that this viscosity is
*       vsa(,,is) = mp*am(is)*na(,,is)*dfvsa
*     If cfvsa(7,is) is not zero, the vsa0 array is flux-scaled with
*     a scaling factor of cfvsa(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfvsa(4:6,) are not used in the default version of the code.
*      cfhci - (0:7,0:ns-1) real*8 array.
*     For 0.le.is.lt.ns, cfhci(,is) are the coefficients in the
*     local expression
*       dfhca(,,is) = cfhci(0,is)+cfhci(1,is)/(ne(,)/1.0e20)+
*         cfhci(2,is)*df0(,,is)+
*         cfhci(3,is)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous atom heat conductivity associated
*     with species (is), so that the atom heat conduction coefficient
*     contains a term dfhca(,,is)*na(,,is).
*     If cfhci(7,is) is not zero, the hci arrays are flux-scaled with
*     a scaling factor of cfhci(7,is)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfhci(4:6,) are not used in the default version of the code.
*      cfhce - (0:7) real*8 array.
*     cfhce() are the coefficients in the local expression
*       dfhce(,) = cfhce(0)+cfhce(1)/(ne(,)/1.0e20)+
*         cfhce(3)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the anomalous electron heat conductivity, so that
*     the electron heat conduction coefficient is dfhce(,)*ne(,).
*     If cfhce(7) is not zero, the hce0 array is flux-scaled with
*     a scaling factor of cfhce(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfhce(2) and (4:6) are not used in the default version of the code.
*      cfsig - (0:7) real*8 array.
*     cfsig() are the coefficients in the local expression
*       dfsig(,) = cfsig(0)+cfsig(1)/(ne(,)/1.0e20)+
*         cfsig(3)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     electrical conductivity, which is
*       sig(,) = qe*ne(,)*dfsig(,).
*     If cfsig(7) is not zero, the sig0 array is flux-scaled with
*     a scaling factor of cfsig(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfsig(2) and (4:6) are not used in the default version of the code.
*      cfalf - (0:7) real*8 array.
*     cfalf() are the coefficients in the local expression
*       dfalf(,) = cfalf(0)+cfalf(1)/(ne(,)/1.0e20)+
*         cfalf(3)*te(,)/ev/abs(bb(,,3))/16
*     that specifies the diffusivity corresponding to the anomalous
*     thermo-electric coefficient, which is
*       alf(,) = ne(,)*sqrt(qe/te(,))*dfalf(,).
*     If cfalf(7) is not zero, the alf0 array is flux-scaled with
*     a scaling factor of cfalf(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfalf(2) and (4:6) are not used in the default version of the code.
*      cflim - (0:7) real*8 array.
*     cflim contains coefficients for flux limiting. cflim(0)=cflme
*     specifies the flux limit factor for poloidal electron heat
*     transport. cflim(1)=cflmi specifies the flux limit factor for
*     poloidal all atom heat transport. cflim(2)=cflmv specifies the
*     flux limit factor for poloidal atom momentum transport.
*     Zero settings turn off these flux limits.
*     cflim(3:7) are not used in the default version of the code.
*
*     Description of some local variables (including arrays in w(0:*):
*
*     lblcp - character*120.
*     lblcp identifies the instance of the program employed to prepare
*     the table of physics parameters.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ahdr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.nout(0).and.1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..test dimensions
      call xertst (1.le.nx, 'faulty argument nx')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

* ..initialise default physics parameters
*   ..obtain lblcp
      call b2xxid ('b2ah', lblcp(1:60))
      call b2rflb (ninp(0), lblcp(61:120))
*   ..read physics parameters
      call b2rfcp (ninp(0), ns)
*   ..test parameters
      call b2xvcp (nx, ns)

* ..output
*   ..write label
      write (nout(0),'(/2x,a)') 'label for default physics parameters:'
      write (nout(0),'(/4x,a/4x,a)') lblcp(1:60), lblcp(61:120)
*   ..test /b2cmpa/, /b2cmpb/, /b2cmpt/
      call b2xvcp (nx, ns)
*   ..write physics parameters
      call b2wfcp (nout(0), 1, ns)
*   ..write version number
      call cfverw (nout(1), newversion)
*   ..save dimensions
      idum(0) = ns
      call cfwuin (nout(1), 1, idum, 'ns')
*   ..save label
      call cfwuch (nout(1), 120, lblcp, 'label')
*   ..save physics common
      call b2wucp (nout(1), newversion, nx, ns)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ahdr

      end subroutine b2ahdr
