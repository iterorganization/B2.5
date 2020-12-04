*-----------------------------------------------------------------------
*.specification

      subroutine b2trql (nx, ny, ns, gs, vol, hx, hy, qz, qc, pbs,
     &  na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_transport_namelist
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  gs(-1:nx,-1:ny,0:2),
     &  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &  hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1),
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), ne2(-1:nx,-1:ny),
     &  qc(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  chvemx(-1:nx,-1:ny), chvimx(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRQL computes the flux limit velocities.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the following flux limit velocities:
*     - chvemx: electron heat transport.
*     - chvimx: ion heat transport.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2
      save ncall
      real (kind=R8) ::
     *   cfl(-1:nx,-1:ny)
      real (kind=R8) ::
     *   nirm(-1:nx,-1:ny)
*   ..procedures
      real (kind=R8) ::
     *   smin, smax
      external xertst, smin, smax
      external b2tlhe, b2tlhi, b2tlmv, b2xpnr, b2xvsg
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The workspace array w is used as a stack for storage of dynamic
*     real*8 arrays. An array 'name' in w(0:*) will be stored starting at
*     location w(l'name').
*
*     Description of some local variables (including arrays in w(0:*):
*
*     nirm - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), nirm(ix,iy) specifies the atom
*     density weighted by 1/sqrt(am) on the (ix,iy) cell. It will hold
*     that 0.lt.nirm(,).
*     The defining expression for nirm is
*       nirm(,) = (sum is : 0.le.is.lt.ns : na(,,is)/sqrt(am(is))).
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2trql')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol, hx, hy
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*    ..test range of qz, qc
       call xertst (-1.0_R8.lt.smin(n2,qz(-1,-1,0),1).and.
     &   smax(n2,qz(-1,-1,0),1).lt.1.0_R8,
     &   'faulty argument range qz(,,0)')
       call xertst (0.0_R8.lt.smin(n2,qz(-1,-1,1),1).and.
     &   smax(n2,qz(-1,-1,1),1).le.1.0_R8,
     &   'faulty argument range qz(,,1)')
       call xertst (0.0_R8.lt.smin(n2,qc,1).and.
     &   smax(n2,qc,1).le.1.0_R8,
     &   'faulty argument range qc')
*    ..test sign of na, te, ti, ne2
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2, ne2, 1, 'ne2', '.gt.')
      endif
*   ..compute nirm (atom density weighted by 1/sqrt(am))
      call b2xpnr (nx, ny, ns, am, na, nirm)

* ..compute flux limit velocities
*   ..compute electron heat flux limit velocity
      call cfl_tanh_profile(nx, ny, cflim(0),
     .  cflme_core, cflme_sol, cflme_tanh_a, cflme_tanh_b, cfl)
      call b2tlhe (nx, ny, me, cfl, gs, pbs(-1,-1,0), qc, ne, te,
     &  chvemx)
*   ..compute all atom heat flux limit velocity
      call cfl_tanh_profile(nx, ny, cflim(1),
     .  cflmi_core, cflmi_sol, cflmi_tanh_a, cflmi_tanh_b, cfl)
      call b2tlhi (nx, ny, mp, cfl, gs, pbs(-1,-1,0), qc, ni,
     &  nirm, ti, chvimx)

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2trql

      end subroutine b2trql
