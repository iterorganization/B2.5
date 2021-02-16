*-----------------------------------------------------------------------
*.specification

      subroutine b2rups (nget, nx, ny, ns,
     &  ne, na, ua, uadia, te, ti, po,
     &  fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc,
     &  b2fstat_version)
      use b2mod_types
      use b2mod_rates
      use b2mod_time
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nget, nx, ny, ns
      character*10 b2fstat_version
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  ne(-1:nx,-1:ny), uadia(-1:nx,-1:ny,0:1,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     1  fch_32(-1:nx,-1:ny,0:1), fch_52(-1:nx,-1:ny,0:1),
     2  fch_p(-1:nx,-1:ny,0:1), kinrgy(-1:nx,-1:ny,0:ns-1),    !srv 10.06.08
     &  fna_mdf(-1:nx,-1:ny,0:1,0:ns-1), fhe_mdf(-1:nx,-1:ny,0:1),
     &  fhi_mdf(-1:nx,-1:ny,0:1), fna_fcor(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_nodrift(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_he(-1:nx,-1:ny,0:1,0:ns-1),
     &  fnaPSch(-1:nx,-1:ny,0:1,0:ns-1), fhePSch(-1:nx,-1:ny,0:1),
     &  fhiPSch(-1:nx,-1:ny,0:1), fna_eir(-1:nx,-1:ny,0:1,0:ns-1),
     &  fne_eir(-1:nx,-1:ny,0:1), fhe_eir(-1:nx,-1:ny,0:1),
     &  fhi_eir(-1:nx,-1:ny,0:1),
     &  fna_32(-1:nx,-1:ny,0:1,0:ns-1), fna_52(-1:nx,-1:ny,0:1,0:ns-1),
     &  fni_32(-1:nx,-1:ny,0:1), fni_52(-1:nx,-1:ny,0:1),
     &  fne_32(-1:nx,-1:ny,0:1), fne_52(-1:nx,-1:ny,0:1),
     &  fchdia(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1),
     &  fchvispar(-1:nx,-1:ny,0:1), fchvisper(-1:nx,-1:ny,0:1),
     &  fchvisq(-1:nx,-1:ny,0:1), fchinert(-1:nx,-1:ny,0:1),
     &  vaecrb(-1:nx,-1:ny,0:1,0:ns-1), vadia(-1:nx,-1:ny,0:1,0:ns-1),
     &  wadia(-1:nx,-1:ny,0:1,0:ns-1),
     &  veecrb(-1:nx,-1:ny,0:1), vedia(-1:nx,-1:ny,0:1),
     &  floe_noc(-1:nx,-1:ny,0:1), floi_noc(-1:nx,-1:ny,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2RUPS reads the state of the plasma from an un*formatted file.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine reads from unit (nget) the state of the plasma:
*     fields na, ua, te, ti, po, and fluxes fna, fhe, fhi, fch.
*     Unit nget must be connected to an un*formatted file. See the
*     code to discern the required format of the input file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nget - integer, input.
*     nget specifies a unit number for un*formatted input.
*
*     nx, ny, ns - integer, input.
*     See the calling routine for description.
*
*     na - (-1:nx,-1:ny,0:ns-1) real*8 array, output.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), na(ix,iy,is) specifies
*     the density of atomic species (is) on the (ix,iy) cell.
*     It will hold that 0.lt.na(,,).
*
*     ne - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), ne(ix,iy) specifies
*     the electron density on the (ix,iy) cell.
*     It will hold that 0.lt.ne(,).
*
*     ua - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), ua(ix,iy,is) specifies
*     the parallel velocity of atomic species (is) on the (ix,iy) cell.
*     The sign of the parallel velocity indicates whether it is flowing
*     in the same direction as the local poloidal magnetic field.
*
*     uadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), uadia(ix,iy,0,is) specifies
*     the diamagnetic drift velocity of species (is) on the (ix,iy) cell, and
*     uadia(ix,iy,1,is) the radial drift velocity.
* 
*     vadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     (Let vadiax(,,)=vadia(,,0,) and vadiay(,)=vadia(,,1,).)
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), vadiax(ix,iy,is) specifies
*     the effective diamagnetic drift velocity of species (is) in cell (ix,iy)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field).
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
*     perpendicular to the magnetic field).
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
*     perpendicular to the magnetic field).
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
*     perpendicular to the magnetic field).
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
*     perpendicular to the magnetic field).
*     For (ix,iy) in (-1:nx,-1:ny), veecrby(ix,iy) specifies
*     the ExB drift velocity of electrons in cell (ix,iy)
*     in the radial direction.
*     veecrb is computed at cell centers if 'b2tfnb_drift_style' is '0'
*     veecrb is computed at cell faces otherwise (default)
*
*     te - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), te(ix,iy) specifies the electron
*     temperature on the (ix,iy) cell. It will hold that 0.lt.te(,).
*
*     ti - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), ti(ix,iy) specifies the all atom
*     temperature on the (ix,iy) cell. It will hold that 0.lt.ti(,).
*
*     po - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), po(ix,iy) specifies the electric
*     potential on the (ix,iy) cell.
*
*     fna - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array, output.
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
*     fhe - (-1:nx,-1:ny,0:1) real*8 array, output.
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
*     fhi - (-1:nx,-1:ny,0:1) real*8 array, output.
*     (Let fhix(,)=fhi(,,0) and fhiy(,)=fhi(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhix(ix,iy) specifies the all atom
*     heat flux through the face between the (ix,iy) cell and its left
*     neighbor. For ix.eq.-1, fhix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhiy(ix,iy) specifies the all atom
*     heat flux through the face between the (ix,iy) cell and its bottom
*     neighbor. For iy.eq.-1, fhiy(-1:nx,iy) will hold 0.
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
*     fch - (-1:nx,-1:ny,0:1) real*8 array, output.
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
*     kinrgy - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), kinrgy(ix,iy,is) specifies
*     the kinetic energy of a particle of species (is) in cell (ix,iy). 
*
*     tim - real*8.
*     The simulation time (in seconds).
*     It will hold that 0.le.tim.
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer n2, is
      real (kind=R8) :: time(1), zero(-1:nx,-1:ny)
*   ..procedures
      external xertst, cfrure
      external b2xvsg, b2xpne

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2rups')
*   ..set n2 and zero
      n2 = (nx+2)*(ny+2)
      zero = 0.0_R8
*   ..test nget
      call xertst (1.le.nget, 'faulty argument nget')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')

* ..read data
*   ..read plasma state
      call cfrure (nget, n2*ns, na, 'na')
      if(b2fstat_version.ge.'01.001.024') then
        call cfrure (nget, n2, ne, 'ne')
      else
        do is = 0, ns-1
          rza(-1:nx,-1:ny,is) = (zamin(is)+zamax(is))/2.0_R8
          rz2(-1:nx,-1:ny,is) = ((zamin(is)+zamax(is))/2.0_R8)**2
        end do
        call b2xpne (nx, ny, ns, rza, na, zero, ne)
      endif
      call cfrure (nget, n2*ns, ua, 'ua')
      if(b2fstat_version.ge.'01.001.000') then                  !xpb
        call cfrure (nget, n2*2*ns, uadia, 'uadia')             !xpb
      else                                                      !dpc
        uadia=0.0_R8                                            !dpc
      endif                                                     !xpb
      call cfrure (nget, n2, te, 'te')
      call cfrure (nget, n2, ti, 'ti')
      call cfrure (nget, n2, po, 'po')
*   ..read fluxes
      call cfrure (nget, n2*2*ns, fna, 'fna')
      call cfrure (nget, n2*2, fhe, 'fhe')
      call cfrure (nget, n2*2, fhi, 'fhi')
      call cfrure (nget, n2*2, fch, 'fch')
      if(b2fstat_version.ge.'01.001.005') then                  !xpb
        call cfrure (nget, n2*2, fch_32, 'fch_32')
        call cfrure (nget, n2*2, fch_52, 'fch_52')
      else                                                      !dpc
        fch_32(:,:,0)=fch(:,:,0)
        fch_32(:,:,1)=0.0_R8
        fch_52(:,:,0)=0.0_R8
        fch_52(:,:,1)=fch(:,:,1)
      endif
      if(b2fstat_version.ge.'01.001.009') then
        call cfrure (nget, n2*ns, kinrgy, 'kinrgy')
      else
        do is = 0, ns-1
          kinrgy(:,:,is)=0.5_R8*am(is)*mp*ua(:,:,is)**2
        enddo
      endif
      call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
*   ..read time
      if(b2fstat_version.ge.'01.001.018') then                  !xpb
        call cfrure (nget, 1, time, 'time')
      else
        time(1)=0.0_R8
      endif
      tim=time(1)
      call xertst (0.0_R8.le.tim, 'faulty input time')
*   ..read parallel current
      if(b2fstat_version.ge.'01.001.069') then                   !srv 09.09.11
        call cfrure (nget, n2, fch_p(:,:,0), 'fch_p')
        fch_p(:,:,1)=0.0_R8                                      !srv 19.03.12
      elseif(b2fstat_version.ge.'01.001.025' .and.
     &       b2fstat_version.le.'01.001.031') then  ! clash between 5.0/1 and 5.2
! cxpb This code looks strange but the idea is that cfrure_opt will fail and rewind
! cxpb if the size of fch_p is n2 and not n2*2
        fch_p(:,:,1)=0.0_R8
        call cfrure_opt (nget, n2*2, fch_p, 'fch_p')
        write(6,*) 'Recovering by reading fch_p with array size ', n2
        call cfrure (nget, n2, fch_p(:,:,0), 'fch_p')
      elseif(b2fstat_version.ge.'01.001.032' .and.
     &       b2fstat_version.le.'01.001.068') then               !xpb srv 10.06.08 {
        call cfrure (nget, n2*2, fch_p, 'fch_p')
      elseif (b2fstat_version.ge.'01.001.022' .and.
     &        b2fstat_version.ge.'01.001.024') then              !srv 10.06.08 }
        call cfrure (nget, n2, fch_p(:,:,0), 'fch_p')
        fch_p(:,:,1)=0.0_R8                                      !srv 19.03.12
      else
        fch_p(:,:,0) = fch(:,:,0)     ! this ensures the old behaviour is kept
        fch_p(:,:,1) = 0.0_R8         ! this ensures the old behaviour is kept
      endif                                                      !srv 10.06.08 }
*   ..adding various fluxes and velocities needed by advanced SOLPS-ITER features
      if(b2fstat_version.ge.'03.000.005') then
        call cfrure (nget, n2*2*ns, fna_mdf, 'fna_mdf')
        call cfrure (nget, n2*2, fhe_mdf, 'fhe_mdf')
        call cfrure (nget, n2*2, fhi_mdf, 'fhi_mdf')
        call cfrure (nget, n2*2*ns, fna_fcor, 'fna_fcor')
        call cfrure (nget, n2*2*ns, fna_nodrift, 'fna_nodrift')
        call cfrure (nget, n2*2*ns, fna_he, 'fna_he')
        call cfrure (nget, n2*2*ns, fnaPSch, 'fnaPSch')
        call cfrure (nget, n2*2, fhePSch, 'fhePSch')
        call cfrure (nget, n2*2, fhiPSch, 'fhiPSch')
        call cfrure (nget, n2*2*ns, fna_eir, 'fna_eir')
        call cfrure (nget, n2*2, fne_eir, 'fne_eir')
        call cfrure (nget, n2*2, fhe_eir, 'fhe_eir')
        call cfrure (nget, n2*2, fhi_eir, 'fhi_eir')
        call cfrure (nget, n2*2*ns, fna_32, 'fna_32')
        call cfrure (nget, n2*2*ns, fna_52, 'fna_52')
        call cfrure (nget, n2*2, fni_32, 'fni_32')
        call cfrure (nget, n2*2, fni_52, 'fni_52')
        call cfrure (nget, n2*2, fne_32, 'fne_32')
        call cfrure (nget, n2*2, fne_52, 'fne_52')
        call cfrure (nget, n2*2, fchdia, 'fchdia')
        call cfrure (nget, n2*2, fchin, 'fchin')
        call cfrure (nget, n2*2, fchvispar, 'fchvispar')
        call cfrure (nget, n2*2, fchvisper, 'fchvisper')
        call cfrure (nget, n2*2, fchvisq, 'fchvisq')
        call cfrure (nget, n2*2, fchinert, 'fchinert')
        call cfrure (nget, n2*2*ns, vaecrb, 'vaecrb')
        call cfrure (nget, n2*2*ns, vadia, 'vadia')
        call cfrure (nget, n2*2*ns, wadia, 'wadia')
        call cfrure (nget, n2*2, veecrb, 'veecrb')
        call cfrure (nget, n2*2, vedia, 'vedia')
        call cfrure (nget, n2*2, floe_noc, 'floe_noc')
        call cfrure (nget, n2*2, floi_noc, 'floi_noc')
      else
        fna_mdf = fna
        fhe_mdf = fhe
        fhi_mdf = fhi
        fna_fcor = fna
        fna_nodrift = fna
        fna_he = fna
        fnaPSch = 0.0_R8
        fhePSch = 0.0_R8
        fhiPSch = 0.0_R8
        fna_eir = fna
        fhe_eir = fhe
        fhi_eir = fhi
        fna_32(:,:,0,:)=fna(:,:,0,:)
        fna_32(:,:,1,:)=0.0_R8
        fna_52(:,:,0,:)=0.0_R8
        fna_52(:,:,1,:)=fna(:,:,1,:)
        call b2xpfi (nx, ny, ns, fna_32, fni_32)
        call b2xpfi (nx, ny, ns, fna_52, fni_52)
        call b2xpfe (nx, ny, ns, 0, qe,
     1     rza, zero, fna_32, zero, fch_32, fne_32)
        call b2xpfe (nx, ny, ns, 0, qe,
     1     rza, zero, fna_52, zero, fch_52, fne_52)
        call b2xpfe (nx, ny, ns, 0, qe,
     1     rza, zero, fna_eir, zero, fch, fne_eir)
        fchdia = 0.0_R8
        fchin = 0.0_R8
        fchvispar = 0.0_R8
        fchvisper = 0.0_R8
        fchvisq = 0.0_R8
        fchinert = 0.0_R8
        vaecrb = 0.0_R8
        vadia = 0.0_R8
        wadia = 0.0_R8
        veecrb = 0.0_R8
        vedia = 0.0_R8
        floe_noc = 0.0_R8
        floi_noc = 0.0_R8
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2rups

      end subroutine b2rups
