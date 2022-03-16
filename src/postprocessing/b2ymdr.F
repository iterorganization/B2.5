*-----------------------------------------------------------------------
*.specification

      subroutine b2ymdr (ninp, nout, nx, ny, ns)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_work
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:1), nout(0:1), nx, ny, ns
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YMDR is a postprocessor for b2. It is meant to produce movie
*     output.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine is a postprocessor for the b2 code. It displays the
*     evolution of the plasma state. This routine does not require any
*     metric information.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:1) integer array, input.
*     ninp specifies the input unit numbers as follows.
*     ninp(0): formatted; provides output control parameters.
*     ninp(1): un*formatted; provides the plasma state evolution.
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
      integer n2, itim, idum(0:9), iftp
      character lblmn*120, idout*80, id*8
      real (kind=R8) ::
     *   tim, dtim, rdum(0:9)
*   ..procedures
      logical streql, cfendf
      external xertst, xerrab, streql, gacwk,
     &  gdawk, cfruch, cfruin, cfrure, cfendf, lbfilx
      external b2wgdt, b2ruzd

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
*     delna, delua, delte, delti, delpo - real*8 array.
*     These arrays describe the change in the plasma state over one
*     timestep.
*
*     lblmn - character*120
*     lblmn identifies the instance of the program employed for the main
*     calculation.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ymdr')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.nout(0).and.
     &             1.le.nout(1), 'faulty argument ninp, nout')
*   ..test dimensions
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')
*   ..read output control
      read (ninp(0),'(a8)',err=91) id
      call xertst (streql(id,'*idout'),
     &  'faulty input instead of *idout, id='//id)
      read (ninp(0),*,err=92) iftp, idout
*   ..read parameters
      call cfruch (ninp(1), 120, lblmn, 'label')
      call b2ruzd (ninp(1), b2fmovie_version,
     &  ns, zamin, zamax, zn, am, .true.)
*   ..write label
      write (nout(0),'(/2x,a)') 'label for this calculation:'
      write (nout(0),'(/4x,a/4x,a)') lblmn(1:60), lblmn(61:120)
*   ..write output control
      write (nout(0),'(/2x,a,2x,i3/2x,a,2x,a)')
     &  'style of contour plots:', iftp,
     &  'idout:', idout

* ..produce (substitute-)movie
*   ..activate graphics on nout(1)
      call gacwk (nout(1))
*   ..mark contour style
      call lbfilx (iftp)
*   ..main loop
    1 if (.not.cfendf(ninp(1))) then
*    ..get itim
       call cfruin (ninp(1), 1, idum, 'itim')
       itim = idum(0)
*    ..get tim, dtim
       call cfrure (ninp(1), 2, rdum, 'tim,dtim')
       tim = rdum(0)
       dtim = rdum(1)
*    ..echo itim, tim, dtim
       write (nout(0),'(2x,a,i5,5x,a,1pe11.4,5x,a,e9.2)')
     &   'itim:', itim, 'tim:', tim, 'dtim:', dtim
*    ..get change in plasma state
       call cfrure (ninp(1), n2*ns, delna, 'delna')
       call cfrure (ninp(1), n2*ns, delua, 'delua')
       call cfrure (ninp(1), n2, delte, 'delte')
       call cfrure (ninp(1), n2, delti, 'delti')
       call cfrure (ninp(1), n2, delpo, 'delpo')
*    ..plot change in state
       call b2wgdt (nout(0), idout, nx, ny, ns,
     &   delne, delni, delna, delua, delte, delti, delpo)
*    ..loop
       goto 1
      endif
*   ..deactivate graphics on nout(1)
      call gdawk (nout(1))

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read heading *idout')
   92 call xerrab ('error trying to read output control')

*-----------------------------------------------------------------------
*.end b2ymdr

      end subroutine b2ymdr
