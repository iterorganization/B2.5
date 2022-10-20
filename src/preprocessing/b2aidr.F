*-----------------------------------------------------------------------
*.specification

      subroutine b2aidr (ninp, nout, nx, ny, ns)
      use b2mod_types
      use b2mod_rates
      use b2mod_geo
      use b2mod_time
      use b2mod_work
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2us_prep
      use b2us_io
      use b2mod_ipmain
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:3), nout(0:1), nx, ny, ns                  !srv 22.12.09
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AIDR prepares and writes to an external file an initial plasma
*     state for B2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine prepares and writes to an external file an initial
*     plasma state for B2.
*     The initial plasma state is characterised by constant densities
*     and temperatures and vanishing flow.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:2) integer array, input.
*     ninp specifies the input unit numbers, as follows.
*     ninp(0): formatted; not used here.
*     ninp(1): formatted; provides data for the initial state.
*     ninp(2): un*formatted; provides the geometry.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides summary printed output.
*     nout(1): un*formatted; provides the initial state.
*
*     nx, ny - integers, input.
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
      integer n2, is, js, i, is0, idum(0:9), ninp_us, nout_us, io
      integer, save :: read_b2fstate_file = 0                       !srv 22.12.09
      real (kind=R8) :: t0, t1                                      !srv 22.12.09
      type(mapping) :: mpg
      type(geometry) :: geo
      type(B2State) :: state
      real (kind=R8), allocatable :: zcdum(:,:),
     &  na1(:,:,:), ua1(:,:,:),  ne1(:,:), uadia1(:,:,:,:,:),
     &  te1(:,:), ti1(:,:), po1(:,:), fna1(:,:,:,:,:), fhe1(:,:,:,:),
     &  fhi1(:,:,:,:), fhn1(:,:,:,:), fch1(:,:,:,:),
     &  fch_321(:,:,:,:), fch_521(:,:,:,:),
     &  fch_p1(:,:,:,:), kinrgy1(:,:,:),
     &  fna_mdf1(:,:,:,:,:), fhe_mdf1(:,:,:,:), fhi_mdf1(:,:,:,:),
     &  fna_fcor1(:,:,:,:,:), fna_nodrift1(:,:,:,:,:),
     &  fna_he1(:,:,:,:,:),
     &  fnaPSch1(:,:,:,:,:), fhePSch1(:,:,:,:), fhiPSch1(:,:,:,:),
     &  fna_eir1(:,:,:,:,:), fne_eir1(:,:,:,:),
     &  fhe_eir1(:,:,:,:), fhi_eir1(:,:,:,:),
     &  fna_321(:,:,:,:,:), fna_521(:,:,:,:,:),
     &  fni_321(:,:,:,:), fni_521(:,:,:,:),
     &  fne_321(:,:,:,:), fne_521(:,:,:,:),
     &  fchdia1(:,:,:,:), fchin1(:,:,:,:),
     &  fchvispar1(:,:,:,:), fchvisper1(:,:,:,:),
     &  fchvisq1(:,:,:,:), fchinert1(:,:,:,:),
     &  vaecrb1(:,:,:,:,:), vadia1(:,:,:,:,:), wadia1(:,:,:,:,:),
     &  veecrb1(:,:,:,:), vedia1(:,:,:,:),
     &  floe_noc1(:,:,:,:), floi_noc1(:,:,:,:),
     &  tn1(:,:), kt1(:,:), zt1(:,:)
      real (kind=R8) :: zpdum(0:ns-1,0:3), zero(-1:nx,-1:ny)
      character lblps*120, lblgm*120, id*8, id0*8, line*500, ligne*500
      character cnamip*80, cvalip*80                                !srv 22.12.09
      real (kind=R8) ::
     *  naini(0:ns-1), ttini, na_min
      logical oldversion, ex
      character*256 filename
*   ..procedures
      logical streql, is_comment
      logical ltst                                                  !srv 22.12.09
      external xertst, xerrab, sfill, streql, cfwuch, cfwuin, is_comment
      external b2xvsg, b2xxid, b2rflb, b2wups, b2rugm, b2xpne, b2wuzd,
     &  b2ruzd, b2rups
      ltst(t0,t1) = abs(t0-t1).le.1.0e-6_R8*(abs(t0)+abs(t1))       !srv 22.12.09
*  ..initialisation
      save na_min
      data na_min/1.0e4_R8/

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
*     See the calling routine for description.
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
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), ua(ix,iy,is) specifies
*     the parallel velocity of atomic species (is) on the (ix,iy) cell.
*     The sign of the parallel velocity indicates whether it is flowing
*     in the same direction as the local poloidal magnetic field.
*
*     uadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), uadia(ix,iy,0,is) specifies
*     the diamagnetic drift velocity of species (is) on the (ix,iy) cell,
*     projected onto the poloidal direction, and
*     uadia(ix,iy,1,is) the radial drift velocity.
*     uadia is computed at cell centers if 'b2tfnb_drift_style' is '0'
*     uadia is computed at cell faces otherwise (default)
*
*     vadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let vadiax(,,)=vadia(,,0,) and vadiay(,)=vadia(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vadiax(ix,iy,is) specifies
*     the effective diamagnetic drift velocity of species (is) in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but
*     projected onto the poloidal direction.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vadiay(ix,iy,is) specifies
*     the effective diamagnetic drift velocity of species (is) in cell (ix,iy)
*     in the radial direction.
*     vadia is computed at cell centers if 'b2tfnb_drift_style' is '0'
*     vadia is computed at cell faces otherwise (default)
*
*     wadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let wadiax(,,)=wadia(,,0,) and wadiay(,)=wadia(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), wadiax(ix,iy,is) specifies
*     the total diamagnetic drift velocity of species (is) in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but
*     projected onto the poloidal direction.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), wadiay(ix,iy,is) specifies
*     the total diamagnetic drift velocity of species (is) in cell (ix,iy)
*     in the radial direction.
*     wadia is computed at cell centers if 'b2tfnb_drift_style' is '0'
*     wadia is computed at cell faces otherwise (default)
*
*     vedia - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let vediax(,)=vedia(,,0) and vediay(,)=vedia(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), vediax(ix,iy) specifies
*     the effective diamagnetic drift velocity of electrons in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field),
*     projected onto the poloidal direction.
*     For (ix,iy) in (-1:nx,-1:ny), vediay(ix,iy) specifies
*     the effective diamagnetic drift velocity of electrons in cell (ix,iy)
*     in the radial direction.
*     vedia is computed at cell centers if 'b2tfnb_drift_style' is '0'
*     vedia is computed at cell faces otherwise (default)
*
*     vaecrb - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let vaecrbx(,,)=vaecrb(,,0,) and vaecrby(,)=vaecrb(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vaecrbx(ix,iy,is) specifies
*     the ExB drift velocity of species (is) in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but
*     projected onto the poloidal direction.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vaecrby(ix,iy,is) specifies
*     the ExB drift velocity of species (is) in cell (ix,iy)
*     in the radial direction.
*     vaecrb is computed at cell centers if 'b2tfnb_drift_style' is '0'
*     vaecrb is computed at cell faces otherwise (default)
*
*     veecrb - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let veecrbx(,)=veecrb(,,0) and veecrby(,)=veecrb(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), veecrbx(ix,iy) specifies
*     the ExB drift velocity of electrons in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but
*     projected onto the poloidal direction.
*     For (ix,iy) in (-1:nx,-1:ny), veecrby(ix,iy) specifies
*     the ExB drift velocity of electrons in cell (ix,iy)
*     in the radial direction.
*     veecrb is computed at cell centers if 'b2tfnb_drift_style' is '0'
*     veecrb is computed at cell faces otherwise (default)
*
*     kinrgy - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), kinrgy(ix,iy,is) specifies
*     the kinetic energy of a particle of species (is) in cell (ix,iy).
*     The velocity used includes the parallel velocity as well
*     as the total diamagnetic and ExB drift velocities in both the radial
*     and diamagnetic directions.
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
*     (ix,iy) cell and its left neighbor.
*     For ix.eq.-1, fnax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fnay(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix,iy) cell and its bottom neighbor.
*     For iy.eq.-1, fnay(-1:nx,iy,0:ns-1) will hold 0.
*
*     fna_mdf - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let fnax_mdf(,,)=fna_mdf(,,0,) and fnay_mdf(,,)=fna_mdf(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), fnax_mdf(ix,iy,is) specifies
*     the modified flux of atoms of species (is) through the face between
*     the (ix,iy) cell and its left neighbor.
*     For ix.eq.-1, fnax_mdf(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fnay_mdf(ix,iy,is) specifies
*     the modified flux of atoms of species (is) through the face between
*     the (ix,iy) cell and its bottom neighbor.
*     For iy.eq.-1, fnay_mdf(-1:nx,iy,0:ns-1) will hold 0.
*     fna_mdf is equal fna for runs without drifts.
*
*     fna_fcor - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let fna_fcorx(,,)=fna_fcor(,,0,) and fna_fcory(,,)=fna_fcor(,,1,).)
*     For (ix,iy,is) in (0:nx,-1:ny,0:ns-1), fna_fcorx(ix,iy,is) specifies
*     the flux of atoms of species (is) for calculating momentum transport
*     through the face between the (ix,iy) cell and its left neighbor.
*     For ix.eq.-1, fna_fcorx(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fna_fcory(ix,iy,is) specifies
*     the flux of atoms of species (is) for calculating momentum transport
*     through the face between the (ix,iy) cell and its bottom neighbor.
*     For iy.eq.-1, fna_fcory(-1:nx,iy,0:ns-1) will hold 0.
*     fna_fcor takes into account Coriolis force.
*     fna_fcor is equal fna for runs without drifts.
*
*     fhe - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhex(,)=fhe(,,0) and fhey(,)=fhe(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhex(ix,iy) specifies the electron
*     heat flux through the face between the (ix,iy) cell and its
*     left neighbor. For ix.eq.-1, fhex(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhey(ix,iy) specifies the electron
*     heat flux through the face between the (ix,iy) cell and its
*     bottom neighbor. For iy.eq.-1, fhey(-1:nx,iy) will hold 0.
*
*     fhe_mdf - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhex_mdf(,)=fhe_mdf(,,0) and fhey_mdf(,)=fhe_mdf(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhex_mdf(ix,iy) specifies the modified
*     electron heat flux through the face between the (ix,iy) cell
*     and its left neighbor. For ix.eq.-1, fhex_mdf(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhey_mdf(ix,iy) specifies the modified
*     electron heat flux through the face between the (ix,iy) cell and
*     its bottom neighbor. For iy.eq.-1, fhey_mdf(-1:nx,iy) will hold 0.
*     fhe_mdf is equal fhe for runs without drifts.
*
*     fhi - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhix(,)=fhi(,,0) and fhiy(,)=fhi(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhix(ix,iy) specifies the all atom
*     heat flux through the face between the (ix,iy) cell and its
*     left neighbor. For ix.eq.-1, fhix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhiy(ix,iy) specifies the all atom
*     heat flux through the face between the (ix,iy) cell and its
*     bottom neighbor. For iy.eq.-1, fhiy(-1:nx,iy) will hold 0.
*
*     fhi_mdf - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fhix_mdf(,)=fhi_mdf(,,0) and fhiy_mdf(,)=fhi_mdf(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhix_mdf(ix,iy) specifies the
*     modified all atom heat flux through the face between the (ix,iy) cell
*     and its left neighbor. For ix.eq.-1, fhix_mdf(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhiy_mdf(ix,iy) specifies the modified
*     all atom heat flux through the face between the (ix,iy) cell and
*     its bottom neighbor. For iy.eq.-1, fhiy_mdf(-1:nx,iy) will hold 0.
*
*     fch - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchx(,)=fch(,,0) and fchy(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchx(ix,iy) specifies the electric
*     current through the face between the (ix,iy) cell and its left
*     neighbor. For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchy(ix,iy) specifies the electric
*     current through the face between the (ix,iy) cell and its bottom
*     neighbor. For iy.eq.-1, fchy(-1:nx,iy) will hold 0.
*
*     fch_p - (-1:nx,-1:ny,0:1) real*8 array.
*     For (ix,iy) in (0:nx,-1:ny), fch_p(ix,iy,0) specifies the product of
*     the parallel electric current and the poloidal magnetic field component
*     sign through the face between the (ix,iy) cell and its left neighbor.
*     For ix.eq.-1, fch_p(ix,-1:ny,0) will hold 0.
*     fch_p(-1:nx,-1:ny,1) = 0
*
*     fchdia - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchdiax(,)=fchdia(,,0) and fchdiay(,)=fchdia(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchdiax(ix,iy) specifies the
*     modified diamagnetic current through the face between
*     the (ix,iy) cell and its left neighbor.
*     For ix.eq.-1, fchdiax(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchdiay(ix,iy) specifies the
*     modified diamagnetic current through the face between
*     the (ix,iy) cell and its bottom neighbor.
*     For iy.eq.-1, fchdiay(-1:nx,iy) will hold 0.
*
*     fchvispar - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchvisparx(,)=fchvispar(,,0) and fchvispary(,)=fchvispar(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchvisparx(ix,iy) specifies the
*     electric current driven by parallel viscosity through the face
*     between the (ix,iy) cell and its left neighbor.
*     For ix.eq.-1, fchvisparx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchvispary(ix,iy) specifies the
*     electric current driven by parallel viscosity through the face
*     between the (ix,iy) cell and its bottom neighbor.
*     For iy.eq.-1, fchvispary(-1:nx,iy) will hold 0.
*
*     fchvisper - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchvisperx(,)=fchvisper(,,0) and fchvispery(,)=fchvisper(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), fchvisperx(ix,iy) = 0.
*     For (ix,iy) in (-1:nx, 0:ny), fchvispery(ix,iy) specifies the
*     electric current connected with contribution from
*     the perpendicular viscosity through the face between the (ix,iy) cell
*     and its bottom neighbor.
*     For iy.eq.-1, fchvispery(-1:nx,iy) will hold 0.
*
*     fchvisq - (-1:nx,-1:ny,0:1) real*8 array.
*     fchvisq specifies the electric current produced by the components
*     of the viscosity tensor which is connected with the heat fluxes.
*     (Let fchvisqx(,)=fchvisq(,,0) and fchvisqy(,)=fchvisq(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchvisqx(ix,iy) specifies the current
*     through the face between the (ix,iy) cell and its left neighbor.
*     For ix.eq.-1, fchvisqx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchvisqy(ix,iy) specifies the current
*     through the face between the (ix,iy) cell and its bottom neighbor.
*     For iy.eq.-1, fchvisqy(-1:nx,iy) will hold 0.
*
*     fchinert - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchinertx(,)=fchinert(,,0) and fchinerty(,)=fchinert(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchinertx(ix,iy) specifies the
*     contributions from inertia and gyroviscosity current through
*     the face between the (ix,iy) cell and its left neighbor.
*     For ix.eq.-1, fchinertx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchinerty(ix,iy) specifies the
*     contributions from inertia and gyroviscosity current through
*     the face between the (ix,iy) cell and its bottom neighbor.
*     For iy.eq.-1, fchinerty(-1:nx,iy) will hold 0.
*
*     fchin - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let fchinx(,)=fchin(,,0) and fchiny(,)=fchin(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchinx(ix,iy) specifies the
*     contribution from ion-neutral friction through the face between
*     the (ix,iy) cell and its left neighbor.
*     For ix.eq.-1, fchinx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchiny(ix,iy) specifies the
*     contribution from ion-neutral friction through the face between
*     the (ix,iy) cell and its bottom neighbor.
*     For iy.eq.-1, fchiny(-1:nx,iy) will hold 0.
*
*     lblps - character*120.
*     lblps identifies the instance of the program employed to prepare
*     the initial state.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2aidr')
*   ..set n2 and zero
      n2 = (nx+2)*(ny+2)
      zero = 0.0_R8
*   ..test ninp, nout
      call xertst (1.le.ninp(1).and.1.le.ninp(2).and.
     &             1.le.nout(0).and.1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..test dimensions
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')

* ..obtain initial plasma state
*   ..obtain lblps
      call b2xxid ('b2ai', lblps(1:60))
      call b2rflb (ninp(1), lblps(61:120))
*   ..read and test species parameters
      read (ninp(1),'(a8)',err=91) id
      call xertst (streql(id,'*specs'),
     &  'faulty input id instead of *specs, id='//id)
*   ..read species data
*    ..Two formats are accepted:
*    ..old version: za, zn, am, za2
*    ..new version: zamin, zamax, zn, am
      oldversion = .true.
      do is = 0, ns-1
       read (ninp(1),*,err=92) id0, (zpdum(is,i),i=0,3)
       read (id0,'(i8)',err=93) is0
       call xertst (is0.eq.is, 'faulty input under *specs')
       oldversion = oldversion .and.
     &   nint(zpdum(is,0)**2).eq.nint(zpdum(is,3)) .and.
     &   zpdum(is,1) .gt. 0.0_R8
      enddo
      if (oldversion) then
       do is = 0, ns-1
        is_neutral(is) = zpdum(is,0).eq.0.0_R8
        zamin(is) = zpdum(is,0)
        zamax(is) = zpdum(is,0)
        zn(is) = zpdum(is,1)
        am(is) = zpdum(is,2)
       enddo
      else
       do is = 0, ns-1
        is_neutral(is) = zpdum(is,0).eq.0.0_R8
        zamin(is) = zpdum(is,0)
        zamax(is) = zpdum(is,1)
        zn(is) = zpdum(is,2)
        am(is) = zpdum(is,3)
       enddo
      endif
      nspecies = 1
      is = 0
      partition_number(is) = 0
      do while (is.lt.ns-1)
        if (lnext(is,is+1)) then
          partition_number(is+1) = partition_number(is)+1
        else
          partition_number(is+1) = 0
          nspecies = nspecies + 1
        endif
        is = is + 1
      end do
      do is = 0, ns-1
        call xertst (0.0_R8.le.zamin(is).and.zamin(is).le.zamax(is),
     .    'faulty input range zamin, zamax')
        call xertst (0.0_R8.lt.zn(is), 'faulty input zn')
        call xertst (zamin(is).eq.nint(zamin(is)),
     .    'faulty input zamin; not integer')
        call xertst (zamax(is).eq.nint(zamax(is)),
     .    'faulty input zamax; not integer')
        call xertst (zn(is).eq.nint(zn(is)),
     .    'faulty input zn; not integer')
        call xertst (0.0_R8.lt.am(is), 'faulty input am')
        if (is_neutral(is)) then
          call xertst (0.0_R8.eq.zamax(is),
     .      'no bundling of neutral species allowed')
        endif
      enddo
*   ..read and test naini
      read (ninp(1),'(a8)',err=91) id
      call xertst (streql(id,'*naini'),
     &  'faulty input id instead of *naini, id='//id)
      read (ninp(1),*,err=92) (naini(is),is=0,ns-1)
      call b2xvsg (ns, naini, 1, 'naini', '.gt.')
*   ..read and test ttini
*     (ttini is provided in units of eV)
      read (ninp(1),'(a8)',err=91) id
      call xertst (streql(id,'*ttini'),
     &  'faulty input id instead of *ttini, id='//id)
      read (ninp(1),*,err=92) ttini
      call xertst (0.lt.ttini, 'faulty input ttini')
*   ..read and echo code internal parameters                      !srv 22.12.09 {
      write (nout(0),'(/2x,a)') 'non-default internal parameters:'
    1 continue
       read (ninp(1),'(a)',end=2,err=94) line
       if (.not.is_comment(line)) then
         ligne = line
         call strip_spaces(ligne)
         if (ligne(1:1).eq.'''') then
           read (line,*) cnamip, cvalip
           write (nout(0),'(4x,a,2x,a)') cnamip, trim(cvalip)
           call ipsetc (cnamip, cvalip)
         endif
       endif
       goto 1
    2 continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'
      call ipgeti ('b2aidr_read_b2fstate', read_b2fstate_file)
      call ipgetr ('b2mndr_na_min', na_min)
      call xertst (na_min.gt.0.0_R8,
     & 'faulty internal parameter na_min')
      do is=0,ns-1
        call xertst (naini(is).ge.na_min,
     &   'Initial densities must be larger than na_min!')
      end do

*   ..initialise time
      tim = 0.0_R8
*   ..obtain version number                                          !xpb
      call cfverr (ninp(2), b2fgmtry_version)                        !xpb
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(2), 2, idum, 'nx,ny')
      call xertst (idum(0).eq.nx.and.idum(1).eq.ny,
     &  'faulty input nx, ny from geometry file')
*   ..read header for geometry
      call cfruch (ninp(2), 120, lblgm, 'label')
      call cfruin (ninp(2), 1, idum, 'isymm')
      isymm = idum(0)
*     (isymm is now stored in the b2cmfs common)                     !xpb 10.08.00
*   ..read geometry
      call b2rugm (ninp(2), nx, ny, isymm, crx, cry, fpsi, ffbz,
     &  bb, vol, hx, hy, hc, ht, qz, qc, qs, qac, qas, gs, pbs,
     &  wbbl, wbbr, wbbv, wbbc, ebc,
     &  cell_width, cell_height, gmap)

      if(read_b2fstate_file.ne.0) then
*  .. copy the data of first n=min(ns0,ns1) species
*     from the plasma state file b2fstate for ns0 species
*     into the initial plasma state file b2fstati for ns1 species
       filename='b2fstate'
       call find_file(filename,ex)
       if (ex) then
         write(*,*) ' Opening '//trim(filename)
         call cfopen (ninp(3),filename,'old','un*formatted')
       else
         call xerrab ('Could not find old b2fstate file')
       end if
*   ..obtain version number
       call cfverr (ninp(3), b2fstate_version)
*   ..obtain nx0, ny0, ns0 from plasma state file
       call cfruin (ninp(3), 3, idum, 'nx,ny,ns')
       write(*,*) ' nx0, ny0, ns0', idum(0:2)
       call xertst (idum(0).eq.nx.and.idum(1).eq.ny,
     &  'faulty input nx, ny from plasma state file')
       call xertst (1.le.idum(2),
     &  'faulty input ns from plasma state file')
*   ..read header for plasma initial state
*     (zamin, zamax, zn, am from ninp(3) will be checked against the data
*     already read in)
       allocate(zcdum(0:idum(2)-1,0:4))
       call cfruch (ninp(3), 120, lblps, 'label')
       call b2ruzd (ninp(3), b2fstate_version, idum(2),
     &   zcdum(0,1), zcdum(0,2), zcdum(0,3), zcdum(0,4), .false.)
       do is = 0, min(ns,idum(2))-1
        call xertst (ltst(zcdum(is,1),zamin(is)),
     &  'faulty data for zamin')
        call xertst (ltst(zcdum(is,2),zamax(is)),
     &  'faulty data for zamax')
        call xertst (ltst(zcdum(is,3),zn(is)), 'faulty data for zn')
        call xertst (ltst(zcdum(is,4),am(is)), 'faulty data for am')
       enddo
       allocate(
     *  na1(-1:nx,-1:ny,0:idum(2)-1), ua1(-1:nx,-1:ny,0:idum(2)-1),
     &  ne1(-1:nx,-1:ny), uadia1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  te1(-1:nx,-1:ny), ti1(-1:nx,-1:ny), po1(-1:nx,-1:ny),
     &  fna1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  fhe1(-1:nx,-1:ny,0:1,0:1),
     &  fhi1(-1:nx,-1:ny,0:1,0:1), fhn1(-1:nx,-1:ny,0:1,0:1), 
     &  fch1(-1:nx,-1:ny,0:1,0:1),
     1  fch_321(-1:nx,-1:ny,0:1,0:1), fch_521(-1:nx,-1:ny,0:1,0:1),
     2  fch_p1(-1:nx,-1:ny,0:1,0:1), kinrgy1(-1:nx,-1:ny,0:idum(2)-1),
     &  fna_mdf1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  fhe_mdf1(-1:nx,-1:ny,0:1,0:1), fhi_mdf1(-1:nx,-1:ny,0:1,0:1),
     &  fna_fcor1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  fna_nodrift1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  fna_he1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  fnaPSch1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  fhePSch1(-1:nx,-1:ny,0:1,0:1), fhiPSch1(-1:nx,-1:ny,0:1,0:1),
     &  fna_eir1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  fne_eir1(-1:nx,-1:ny,0:1,0:1), fhe_eir1(-1:nx,-1:ny,0:1,0:1),
     &  fhi_eir1(-1:nx,-1:ny,0:1,0:1),
     &  fna_321(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  fna_521(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  fni_321(-1:nx,-1:ny,0:1,0:1), fni_521(-1:nx,-1:ny,0:1,0:1),
     &  fne_321(-1:nx,-1:ny,0:1,0:1), fne_521(-1:nx,-1:ny,0:1,0:1),
     &  fchdia1(-1:nx,-1:ny,0:1,0:1), fchin1(-1:nx,-1:ny,0:1,0:1),
     &  fchvispar1(-1:nx,-1:ny,0:1,0:1),
     &  fchvisper1(-1:nx,-1:ny,0:1,0:1),
     &  fchvisq1(-1:nx,-1:ny,0:1,0:1), fchinert1(-1:nx,-1:ny,0:1,0:1),
     &  vaecrb1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  vadia1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  wadia1(-1:nx,-1:ny,0:1,0:1,0:idum(2)-1),
     &  veecrb1(-1:nx,-1:ny,0:1,0:1), vedia1(-1:nx,-1:ny,0:1,0:1),
     &  floe_noc1(-1:nx,-1:ny,0:1,0:1), floi_noc1(-1:nx,-1:ny,0:1,0:1),
     &  tn1(-1:nx,-1:ny), kt1(-1:nx,-1:ny), zt1(-1:nx,-1:ny))
*   ..read initial state from the file
       call b2rups (ninp(3), idum(0), idum(1), idum(2),
     &  ne1, na1, ua1, uadia1, te1, ti1, tn1, po1, kt1, zt1,
     &  fna1, fhe1, fhi1, fhn1, fch1, fch_321, fch_521, fch_p1, kinrgy1,
     &  fna_mdf1, fhe_mdf1, fhi_mdf1, fna_fcor1, fna_nodrift1, fna_he1,
     &  fnaPSch1, fhePSch1, fhiPSch1,
     &  fna_eir1, fne_eir1, fhe_eir1, fhi_eir1,
     &  fna_321, fna_521, fni_321, fni_521, fne_321, fne_521,
     &  fchdia1, fchin1, fchvispar1, fchvisper1, fchvisq1, fchinert1,
     &  vaecrb1, vadia1, wadia1, veecrb1, vedia1, floe_noc1, floi_noc1,
     &  b2fstate_version)
       close(ninp(3))
       write(*,*) ' Closed b2fstate'
*   ..test plasma state and fluxes
       call b2xvps_st (idum(0), idum(1), idum(2),                     !srv 20.12.18
     &  ne1, na1, ua1, te1, ti1, po1, fna1, fhe1, fhi1, fch1)
       call b2xvsg (n2*idum(2), kinrgy1, 1, 'kinrgy', '.ge.')
*   ..initialise na, ua, kinrgy, uadia, ne
       do is = 0, min(ns,idum(2))-1
         na(:,:,is)       = na1(:,:,is)
         ua(:,:,is)       = ua1(:,:,is)
         kinrgy(:,:,is)   = kinrgy1(:,:,is)
         uadia(:,:,:,:,is)  = uadia1(:,:,:,:,is)
         vadia(:,:,:,:,is)  = vadia1(:,:,:,:,is)
         wadia(:,:,:,:,is)  = wadia1(:,:,:,:,is)
         vaecrb(:,:,:,:,is) = vaecrb1(:,:,:,:,is)
       enddo
       do is = min(ns,idum(2)), ns-1
        call sfill (n2, naini(is), na(-1,-1,is), 1)
        call sfill (n2, 0.0_R8, ua(-1,-1,is), 1)
        call sfill (n2, 0.0_R8, kinrgy(-1,-1,is), 1)
        call sfill (n2*4, 0.0_R8, uadia(-1,-1,0,0,is), 1)
        call sfill (n2*4, 0.0_R8, vadia(-1,-1,0,0,is), 1)
        call sfill (n2*4, 0.0_R8, wadia(-1,-1,0,0,is), 1)
        call sfill (n2*4, 0.0_R8, vaecrb(-1,-1,0,0,is), 1)
       enddo
       do is = 0, ns-1
        call sfill (n2, (zamin(is)+zamax(is))/2.0_R8, rza(-1,-1,is), 1)
       enddo
       call b2xpne_st (nx, ny, ns, rza, na, zero, ne)                     !srv 20.12.18
*   ..initialise te, ti, tn, po, kt, zt
       te = te1
       ti = ti1
       tn = tn1
       po = po1
       kt = kt1
       zt = zt1
*   ..initialise fluxes
       do is = 0, min(ns,idum(2))-1
         fna(:,:,:,:,is)         = fna1(:,:,:,:,is)
         fna_mdf(:,:,:,:,is)     = fna_mdf1(:,:,:,:,is)
         fna_fcor(:,:,:,:,is)    = fna_fcor1(:,:,:,:,is)
         fna_nodrift(:,:,:,:,is) = fna_nodrift1(:,:,:,:,is)
         fna_he(:,:,:,:,is)      = fna_he(:,:,:,:,is)
         fnaPSch(:,:,:,:,is)     = fnaPSch1(:,:,:,:,is)
         fna_eir(:,:,:,:,is)     = fna_eir1(:,:,:,:,is)
         fna_32(:,:,:,:,is)      = fna_321(:,:,:,:,is)
         fna_52(:,:,:,:,is)      = fna_521(:,:,:,:,is)
       enddo
       do is = min(ns,idum(2)), ns-1
         call sfill (n2*4, 0.0_R8, fna(-1,-1,0,0,is), 1)
         call sfill (n2*4, 0.0_R8, fna_mdf(-1,-1,0,0,is), 1)
         call sfill (n2*4, 0.0_R8, fna_fcor(-1,-1,0,0,is), 1)
         call sfill (n2*4, 0.0_R8, fna_nodrift(-1,-1,0,0,is), 1)
         call sfill (n2*4, 0.0_R8, fna_he(-1,-1,0,0,is), 1)
         call sfill (n2*4, 0.0_R8, fnaPSch(-1,-1,0,0,is), 1)
         call sfill (n2*4, 0.0_R8, fna_eir(-1,-1,0,0,is), 1)
         call sfill (n2*4, 0.0_R8, fna_32(-1,-1,0,0,is), 1)
         call sfill (n2*4, 0.0_R8, fna_52(-1,-1,0,0,is), 1)
       enddo
       vedia     = vedia1
       veecrb    = veecrb1
       fne_eir   = fne_eir1
       fne_32    = fne_321
       fne_52    = fne_521
       fni_32    = fni_321
       fni_52    = fni_521
       fhe       = fhe1
       fhe_mdf   = fhe_mdf1
       fhePSch   = fhePSch1
       fhe_eir   = fhe_eir1
       fhi       = fhi1
       fhn       = fhn1
       fhi_mdf   = fhi_mdf1
       fhiPSch   = fhiPSch1
       fhi_eir   = fhi_eir1
       fch       = fch1
       fch_32    = fch_321
       fch_52    = fch_521
       fch_p     = fch_p1
       fchdia    = fchdia1
       fchin     = fchin1
       fchvispar = fchvispar1
       fchvisper = fchvisper1
       fchvisq   = fchvisq1
       fchinert  = fchinert1
       floe_noc  = floe_noc1
       floi_noc  = floi_noc1
       deallocate( zcdum, na1, ua1, ne1, uadia1, te1, ti1, po1,
     &    fna1, fhe1, fhi1, fhn1, fch1, fch_321, fch_521, fch_p1,
     &    kinrgy1, fna_mdf1, fhe_mdf1, fhi_mdf1,
     &    fna_fcor1, fna_nodrift1, fna_he1,
     &    fnaPSch1, fhePSch1, fhiPSch1,
     &    fna_eir1, fne_eir1, fhe_eir1, fhi_eir1,
     &    fna_321, fna_521, fni_321, fni_521, fne_321, fne_521,
     &    fchdia1, fchin1, fchvispar1, fchvisper1, fchvisq1, fchinert1,
     &    vaecrb1, vadia1, wadia1, veecrb1, vedia1,
     &    floe_noc1, floi_noc1)
      else                                                           !srv 22.12.09 }
*   ..initialise na, ua, kinrgy, uadia, ne
       do is = 0, ns-1
        call sfill (n2, naini(is), na(-1,-1,is), 1)
        call sfill (n2, 0.0_R8, ua(-1,-1,is), 1)
        call sfill (n2, 0.0_R8, kinrgy(-1,-1,is), 1)
        call sfill (n2*4, 0.0_R8, uadia(-1,-1,0,0,is), 1)
        call sfill (n2*4, 0.0_R8, vadia(-1,-1,0,0,is), 1)
        call sfill (n2*4, 0.0_R8, wadia(-1,-1,0,0,is), 1)
        call sfill (n2*4, 0.0_R8, vaecrb(-1,-1,0,0,is), 1)
        call sfill (n2, (zamin(is)+zamax(is))/2.0_R8, rza(-1,-1,is), 1)
       enddo
       call b2xpne_st (nx, ny, ns, rza, na, zero, ne)             !srv 20.12.18
       call sfill (n2*4, 0.0_R8, veecrb, 1)                        !srv 20.12.18
       call sfill (n2*4, 0.0_R8, vedia, 1)                         !srv 20.12.18
*   ..initialise te, ti, tn, po, kt, zt
       call sfill (n2, ttini*ev, te, 1)
       call sfill (n2, ttini*ev, ti, 1)
       call sfill (n2, ttini*ev, tn, 1)
       call sfill (n2, 3.1_R8*ttini, po, 1)
       call sfill (n2, 0.0_R8, kt, 1)
       call sfill (n2, 0.0_R8, zt, 1)
*   ..initialise fluxes
       call sfill (n2*4*ns, 0.0_R8, fna, 1)
       call sfill (n2*4*ns, 0.0_R8, fna_mdf, 1)
       call sfill (n2*4*ns, 0.0_R8, fna_fcor, 1)
       call sfill (n2*4*ns, 0.0_R8, fna_nodrift, 1)
       call sfill (n2*4*ns, 0.0_R8, fna_he, 1)
       call sfill (n2*4*ns, 0.0_R8, fnaPSch, 1)
       call sfill (n2*4*ns, 0.0_R8, fna_eir, 1)
       call sfill (n2*4*ns, 0.0_R8, fna_32, 1)
       call sfill (n2*4*ns, 0.0_R8, fna_52, 1)
       call sfill (n2*4, 0.0_R8, fne_eir, 1)
       call sfill (n2*4, 0.0_R8, fhe_mdf, 1)
       call sfill (n2*4, 0.0_R8, fhe_eir, 1)
       call sfill (n2*4, 0.0_R8, fhePSch, 1)
       call sfill (n2*4, 0.0_R8, fhi_mdf, 1)
       call sfill (n2*4, 0.0_R8, fhi_eir, 1)
       call sfill (n2*4, 0.0_R8, fhiPSch, 1)
       call sfill (n2*4, 0.0_R8, fhe, 1)
       call sfill (n2*4, 0.0_R8, fhi, 1)
       call sfill (n2*4, 0.0_R8, fhn, 1)
       call sfill (n2*4, 0.0_R8, fch, 1)
       call sfill (n2*4, 0.0_R8, fne_32, 1)
       call sfill (n2*4, 0.0_R8, fne_52, 1)
       call sfill (n2*4, 0.0_R8, fni_32, 1)
       call sfill (n2*4, 0.0_R8, fni_52, 1)
       call sfill (n2*4, 0.0_R8, fch_32, 1)
       call sfill (n2*4, 0.0_R8, fch_52, 1)
       call sfill (n2*4, 0.0_R8, fch_p, 1)                           !srv 10.06.08
       call sfill (n2*4, 0.0_R8, fchdia, 1)
       call sfill (n2*4, 0.0_R8, fchin, 1)
       call sfill (n2*4, 0.0_R8, fchvispar, 1)
       call sfill (n2*4, 0.0_R8, fchvisper, 1)
       call sfill (n2*4, 0.0_R8, fchvisq, 1)
       call sfill (n2*4, 0.0_R8, fchinert, 1)
       call sfill (n2*4, 0.0_R8, floe_noc, 1)
       call sfill (n2*4, 0.0_R8, floi_noc, 1)
      endif                                                          !srv 22.12.09
* ..output
*   ..write label
      write (nout(0),'(/2x,a)') 'label for the initial state:'
      write (nout(0),'(/4x,a/4x,a)') lblps(1:60), lblps(61:120)
*   ..write mesh dimensions
      write (nout(0),101) nx, ny
*   ..write species data
      write (nout(0),102) ns,
     &  (is,zamin(is),zamax(is),zn(is),am(is),is=0,ns-1)
*   ..write data on initial state
      write (nout(0),'(/2x,a)') 'Parameters for uniform, stationary'//
     &  ' initial state:'
      do is = 0, ns-1
       write (nout(0),'(2x,a3,i3,a2,2x,1pe11.2)')
     &   'na(',is,'):', naini(is)
      enddo
      write (nout(0),'(7x,a3,2x,1pe11.2)') 'te:', ttini
      write (nout(0),'(7x,a3,2x,1pe11.2)') 'ti:', ttini
      write (nout(0),'(7x,a3,2x,1pe11.2)') 'tn:', ttini
      write (nout(0),'(7x,a3,2x,1pe11.2)') 'po:', 3.1_R8*ttini
*   ..test plasma state
      call b2xvps_st (nx, ny, ns,                                             !srv 20.12.18
     .  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
      call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
*   ..save dimensions
      idum(0) = nx
      idum(1) = ny
      idum(2) = ns
      call cfverw (nout(1), newversion)                                 !xpb
      call cfwuin (nout(1), 3, idum, 'nx,ny,ns')
*   ..save label
      call cfwuch (nout(1), 120, lblps, 'label')
*   ..save /b2cmpa/
      call b2wuzd (nout(1), newversion, ns, zamin, zamax, zn, am)
*   ..save state
      call b2wups (nout(1), nx, ny, ns, ne, na, ua, uadia, te, ti, tn,          !xpb
     &  po, kt, zt, fna, fhe, fhi, fhn, fch, fch_32, fch_52, fch_p,
     &  kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc)

*   ..read unstructured geometry
      open (newunit=ninp_us,file='b2fgmtry_us',iostat=io)
      call cfopen (ninp_us,'b2fgmtry_us','old','un*formatted')
      call cfverr (ninp_us, b2fgmtry_version)
      call cfruin (ninp_us,7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
      mpg%nCi = idum(0) 
      mpg%nCg = idum(1) 
      mpg%nCv = idum(2)
      mpg%nFc = idum(3)
      mpg%nVx = idum(4)
      mpg%nFs = idum(5)
      mpg%nFt = idum(6)
      call cfruin (ninp_us,5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
      mpg%nCmxVx = idum(0)
      mpg%nCmxFc = idum(1)
      mpg%nVmxCv = idum(2)
      mpg%nVmxFc = idum(3)
      mpg%nCmxNv = idum(4)
      call cfruin (ninp_us,1,idum,'isClassicalGrid')
      mpg%isClassicalGrid = idum(0)
      call cfruin (ninp_us,3,idum,'nx,ny,nncut')
      mpg%nx    = idum(0)
      mpg%ny    = idum(1)
      mpg%nncut = idum(2)
      call cfruch (ninp_us, 120, lblgm, 'label')
      call cfruin (ninp_us, 1, idum, 'isymm')
      call alloc_mapping (mpg)
      call alloc_geometry (geo, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs, 
     .  mpg%nFt, mpg%nVmxCv)
      call read_b2fgmtry (ninp_us, mpg, geo)

*   ..convert plasma to unstructured format
      call prep_state (nx, ny, ns,
     .                    mpg%nCv, mpg%nCi, mpg%nFc, mpg, state)
      write (*,*) 'State created'
      call b2xvps (mpg%nCv, mpg%nFc, ns, state%pl, state%dv)

*   ..write unstructured b2fstati
      open (newunit=nout_us,file='b2fstati_us',iostat=io)
      idum(0) = mpg%nCv
      idum(1) = mpg%nFc
      idum(2) = ns
      call cfverw (nout_us, newversion)
      call cfwuin (nout_us, 3, idum, 'nCv,nFc,ns')
      call cfwuch (nout_us, 120, lblps, 'label')
      call b2wuzd (nout_us, newversion, ns, zamin, zamax, zn, am)
      call write_b2fstate (nout_us, mpg%nCv, mpg%nFc, ns, state)

*   ..clean up
      close (ninp_us)
      close (nout_us)
      call dealloc_geometry(geo)
      call dealloc_mapping(mpg)
      call destroyB2State(state)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read heading')
   92 call xerrab ('error trying to read data; id='//id)
   93 call xerrab ('error trying to read species index after '//id)
   94 call xerrab
     1 ('error trying to read non-default internal parameters')      !srv 22.12.09

*-----------------------------------------------------------------------
*.formats

  101 format (/'  The grid contains nx*ny interior cells, ',
     &  'where nx =',i4,' and ny =',i4,'.'/
     &  '  With boundary cells included the size is (nx+2)*(ny+2).')
  102 format (/'  There are ns species, ns =', i3,
     &  '.  The species properties (zamin,zamax,zn,am) are:'/
     &  ('  species ',i3,',  atomic charge range from ',f5.2,
     &   ' to ',f5.2,
     &   ',  nuclear charge ',f5.2,',  atomic mass ',f6.2))

*-----------------------------------------------------------------------
*.end b2aidr

      end subroutine b2aidr
