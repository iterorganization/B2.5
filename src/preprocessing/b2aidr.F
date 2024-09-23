*-----------------------------------------------------------------------
*.specification

      subroutine b2aidr (ninp, nout, nCv, nFc, ns, geo, mpg, st)
      use b2mod_types
      use b2mod_version
      use b2mod_math
      use b2mod_time
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
      integer ninp(0:2), nout(0:1), nCv, nFc, ns                  !srv 22.12.09
*   ..input/output arguments
      type(mapping), intent (inout)  :: mpg
      type(geometry), intent (inout) :: geo
      type(B2State), intent (inout)  :: st
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
*     ninp(0): formatted; provides data for the initial state.
*     ninp(1): un*formatted; provides the geometry.
*     ninp(2): formatted; provides an existing state.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides summary printed output.
*     nout(1): un*formatted; provides the initial state.
*
*     nCv, nFc - integer.
*     nCv and nFc specify the number of control volumes and cell faces
*     in the unstructured mesh.
*     It will hold that 0.lt.nCv and 0.lt.nFc.
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
      integer is, i, is0, ns0, idum(0:9)
      integer, save :: read_b2fstate_file = 0                       !srv 22.12.09
      type(B2State) :: st1
      real (kind=R8), allocatable :: zcdum(:,:)
      real (kind=R8) :: zpdum(0:ns-1,0:3), zero(nCv)
      character lblps*120, lblgm*120, id*8, id0*8, line*500, ligne*500
      character cnamip*80, cvalip*80                                !srv 22.12.09
      real (kind=R8) ::
     *  naini(0:ns-1), ttini, na_min
      logical oldversion, ex
      character*256 filename
*   ..procedures
      logical streql, is_comment
      external xertst, xerrab, sfill, streql, cfwuch, cfwuin, is_comment
      external b2xvsg, b2xxid, b2rflb, b2wups, b2rugm, b2xpne, b2wuzd,
     &  b2ruzd, b2rups
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
*     na - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), na(iCv,is) specifies
*     the density of atomic species (is) on the (iCv) cell.
*     It will hold that 0.lt.na(,).
*
*     ne - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), ne(iCv) specifies
*     the electron density on the (iCv) cell.
*     It will hold that 0.lt.ne().
*
*     ua - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), ua(iCv,is) specifies
*     the parallel velocity of atomic species (is) on the (iCv) cell.
*     The sign of the parallel velocity indicates whether it is flowing
*     in the same direction as the local poloidal magnetic field.
*
*     uadia - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:ns-1), uadia(iFc,0,is) specifies
*     the diamagnetic drift velocity of species (is) on the (iFc) face,
*     projected onto the poloidal direction, and
*     uadia(iFc,1,is) the radial drift velocity.
*
*     vadia - (1:nFc,0:1,0:ns-1) real*8 array.
*     For (iFc,is) in (1:nFc,0:ns-1), vadia(iFc,0,is) specifies
*     the effective diamagnetic drift velocity of species (is) on face (iFc)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but
*     projected onto the poloidal direction.
*     For (iFc,is) in (1:nFc,0:ns-1), vadia(iFc,1,is) specifies
*     the effective diamagnetic drift velocity of species (is) on face (iFc)
*     in the radial direction.
*
*     wadia - (1:nFc,0:1,0:ns-1) real*8 array.
*     (Let wadiax(,)=wadia(,0,) and wadiay(,)=wadia(,1,).)
*     For (iFc,is) in (1:nFc,0:ns-1), wadiax(iFc,is) specifies
*     the total diamagnetic drift velocity of species (is) on face (iFc)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but
*     projected onto the poloidal direction.
*     For (iFc,is) in (1:nFc,0:ns-1), wadiay(iFc,is) specifies
*     the total diamagnetic drift velocity of species (is) on face (iFc)
*     in the radial direction.
*
*     vedia - (1:nFc,0:1) real*8 array.
*     (Let vediax()=vedia(,0) and vediay()=vedia(,1).)
*     For (iFc) in (1:nFc), vediax(iFc) specifies
*     the effective diamagnetic drift velocity of electrons on face (iFc)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field),
*     projected onto the poloidal direction.
*     For (iFc) in (1:nFc), vediay(iFc) specifies
*     the effective diamagnetic drift velocity of electrons on face (iFc)
*     in the radial direction.
*
*     vaecrb - (1:nFc,0:1,0:ns-1) real*8 array.
*     (Let vaecrbx(,)=vaecrb(,0,) and vaecrby(,)=vaecrb(,1,).)
*     For (iFc,is) in (1:nFc,0:ns-1), vaecrbx(iFc,is) specifies
*     the ExB drift velocity of species (is) on face (iFc)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but
*     projected onto the poloidal direction.
*     For (iFc,is) in (1:nFc,0:ns-1), vaecrby(iFc,is) specifies
*     the ExB drift velocity of species (is) on face (iFc)
*     in the radial direction.
*
*     veecrb - (1:nFc,0:1) real*8 array.
*     (Let veecrbx()=veecrb(,0) and veecrby()=veecrb(,1).)
*     For (iFc) in (1:nFc), veecrbx(iFc) specifies
*     the ExB drift velocity of electrons on face (iFc)
*     in the diamagnetic direction (i.e. within the flux surface but
*     perpendicular to the magnetic field), but
*     projected onto the poloidal direction.
*     For (iFc) in (1:nFc), veecrby(iFc) specifies
*     the ExB drift velocity of electrons on face (iFc)
*     in the radial direction.
*
*     kinrgy - (1:nCv,0:ns-1) real*8 array.
*     For (iCv,is) in (1:nCv,0:ns-1), kinrgy(iCv,is) specifies
*     the kinetic energy of a particle of species (is) in cell (iCv).
*     The velocity used includes the parallel velocity as well
*     as the total diamagnetic and ExB drift velocities in both the radial
*     and diamagnetic directions.
*
*     te - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), te(iCv) specifies the electron
*     temperature on the (iCv) cell. It will hold that 0.lt.te().
*
*     ti - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), ti(iCv) specifies the all atom
*     temperature on the (iCv) cell. It will hold that 0.lt.ti().
*
*     po - (1:nCv) real*8 array.
*     For (iCv) in (1:nCv), po(iCv) specifies the electric
*     potential on the (iCv) cell.
*
*     fna - (1:nFc,0:1,0:ns-1) real*8 array.
*     (Let fnax(,)=fna(,0,) and fnay(,)=fna(,1,).)
*     For (iFc,is) in (1:nFc,0:ns-1), fnax(iFc,is) specifies
*     the poloidal flux of atoms of species (is) through the face (iFc).
*     For (iFc,is) in (1:nFc,0:ns-1), fnay(iFc,is) specifies
*     the radial flux of atoms of species (is) through the face (iFc).
*
*     fna_mdf - (1:nFc,0:1,0:ns-1) real*8 array.
*     (Let fnax_mdf(,)=fna_mdf(,0,) and fnay_mdf(,)=fna_mdf(,1,).)
*     For (iFc,is) in (1:nFc,0:ns-1), fnax_mdf(iFc,is) specifies
*     the modified poloidal flux of atoms of species (is) through the face (iFc).
*     For (iFc,is) in (1:nFc,0:ns-1), fnay_mdf(iFc,is) specifies
*     the modified radial flux of atoms of species (is) through the face (iFc).
*     fna_mdf is equal fna for runs without drifts.
*
*     fna_fcor - (1:nFc,0:1,0:ns-1) real*8 array.
*     (Let fna_fcorx(,)=fna_fcor(,0,) and fna_fcory(,)=fna_fcor(,1,).)
*     For (iFc,is) in (1:nFc,0:ns-1), fna_fcorx(iFc,is) specifies
*     the poloidal flux of atoms of species (is) for calculating momentum transport
*     through the face (iFc).
*     For (iFc,is) in (1:nFc,0:ns-1), fna_fcory(iFc,is) specifies
*     the radial flux of atoms of species (is) for calculating momentum transport
*     through the face (iFc).
*     fna_fcor takes into account Coriolis force.
*     fna_fcor is equal fna for runs without drifts.
*
*     fhe - (1:nFc,0:1) real*8 array.
*     (Let fhex()=fhe(,0) and fhey(,)=fhe(,1).)
*     For (iFc) in (1:nFc), fhex(iFc) specifies the poloidal electron
*     heat flux through the face (iFc).
*     For (iFc) in (1:nFc), fhey(iFc) specifies the radial electron
*     heat flux through the face (iFc).
*
*     fhe_mdf - (1:nFc,0:1) real*8 array.
*     (Let fhex_mdf()=fhe_mdf(,0) and fhey_mdf()=fhe_mdf(,1).)
*     For (iFc) in (1:nFc), fhex_mdf(iFc) specifies the modified poloidal
*     electron heat flux through the face (iFc).
*     For (iFc) in (1:nFc), fhey_mdf(iFc) specifies the modified radial
*     electron heat flux through the face (iFc).
*     fhe_mdf is equal fhe for runs without drifts.
*
*     fhi - (1:nFc,0:1) real*8 array.
*     (Let fhix()=fhi(,0) and fhiy()=fhi(,1).)
*     For (iFc) in (1:nFc), fhix(iFc) specifies the poloidal all atom
*     heat flux through the face (iFc).
*     For (iFc) in (1:nFc), fhiy(iFc) specifies the radial all atom
*     heat flux through the face (iFc).
*
*     fhi_mdf - (1:nFc,0:1) real*8 array.
*     (Let fhix_mdf()=fhi_mdf(,0) and fhiy_mdf()=fhi_mdf(,1).)
*     For (iFc) in (1:nFc), fhix_mdf(iFc) specifies the
*     modified poloidal all atom heat flux through the face (iFc).
*     For (iFc) in (1:nFc), fhiy_mdf(iFc) specifies the modified
*     radail all atom heat flux through the face (iFc).
*
*     fch - (1:nFc,0:1) real*8 array.
*     (Let fchx()=fch(,0) and fchy()=fch(,1).)
*     For (iFc) in (1:nFc), fchx(iFc) specifies the poloidal electric
*     current through the face (iFc).
*     For (iFc) in (1:nFc), fchy(iFc) specifies the radial electric
*     current through the face (iFc).
*
*     fch_p - (1:nFc,0:1) real*8 array.
*     For (iFc) in (1:nFc), fch_p(iFc,0) specifies the product of
*     the parallel electric current and the poloidal magnetic field component
*     sign through the face (iFc).
*     fch_p(1:nFc,1) = 0
*
*     fchdia - (1:nFc,0:1) real*8 array.
*     (Let fchdiax()=fchdia(,0) and fchdiay()=fchdia(,1).)
*     For (iFc) in (1:nFc), fchdiax(iFc) specifies the
*     modified poloidal diamagnetic current through the face (iFc).
*     For (iFc) in (1:nFc), fchdiay(iFc) specifies the
*     modified radial diamagnetic current through the face (iFc).
*
*     fchvispar - (1:nFc,0:1) real*8 array.
*     (Let fchvisparx()=fchvispar(,0) and fchvispary()=fchvispar(,1).)
*     For (iFc) in (1:nFc), fchvisparx(iFc) specifies the poloidal
*     electric current driven by parallel viscosity through the face (iFc).
*     For (iFc) in (1:nFc), fchvispary(iFc) specifies the radial
*     electric current driven by parallel viscosity through the face (iFc).
*
*     fchvisper - (1:nFc,0:1) real*8 array.
*     (Let fchvisperx()=fchvisper(,0) and fchvispery()=fchvisper(,1).)
*     For (iFc) in (1:nFc), fchvisperx(iFc) = 0.
*     For (iFc) in (1:nFc), fchvispery(iFc) specifies the radial
*     electric current connected with contribution from
*     the perpendicular viscosity through the face (iFc).
*
*     fchvisq - (1:nFc,0:1) real*8 array.
*     fchvisq specifies the electric current produced by the components
*     of the viscosity tensor which is connected with the heat fluxes.
*     (Let fchvisqx()=fchvisq(,0) and fchvisqy()=fchvisq(,1).)
*     For (iFc) in (1:nFc), fchvisqx(iFc) specifies the poloidal current
*     through the face (iFc).
*     For (iFc) in (1:nFc), fchvisqy(iFc) specifies the radial current
*     through the face (iFc).
*
*     fchinert - (1:nFc,0:1) real*8 array.
*     (Let fchinertx()=fchinert(,0) and fchinerty()=fchinert(,1).)
*     For (iFc) in (1:nFc), fchinertx(iFc) specifies the poloidal
*     contributions from inertia and gyroviscosity current through
*     the face (iFc).
*     For (iFc) in (1:nFc), fchinerty(iFc) specifies the radial
*     contributions from inertia and gyroviscosity current through
*     the face (iFc).
*
*     fchin - (1:nFc,0:1) real*8 array.
*     (Let fchinx()=fchin(,0) and fchiny()=fchin(,1).)
*     For (iFc) in (1:nFc), fchinx(iFc) specifies the poloidal
*     contribution from ion-neutral friction through the face (iFc).
*     For (iFc) in (1:nFc), fchiny(iFc) specifies the radial
*     contribution from ion-neutral friction through the face (iFc).
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
*   ..set zero
      zero = 0.0_R8
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.
     &             1.le.nout(0).and.1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..test dimensions
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')

* ..obtain initial plasma state
*   ..obtain lblps
      call b2xxid ('b2ai', lblps(1:60))
      call b2rflb (ninp(0), lblps(61:120))
*   ..read and test species parameters
      read (ninp(0),'(a8)',err=91) id
      call xertst (streql(id,'*specs'),
     &  'faulty input id instead of *specs, id='//id)
*   ..read species data
*    ..Two formats are accepted:
*    ..old version: za, zn, am, za2
*    ..new version: zamin, zamax, zn, am
      oldversion = .true.
      do is = 0, ns-1
       read (ninp(0),*,err=92) id0, (zpdum(is,i),i=0,3)
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
      read (ninp(0),'(a8)',err=91) id
      call xertst (streql(id,'*naini'),
     &  'faulty input id instead of *naini, id='//id)
      read (ninp(0),*,err=92) (naini(is),is=0,ns-1)
      call b2xvsg (ns, naini, 1, 'naini', '.gt.')
*   ..read and test ttini
*     (ttini is provided in units of eV)
      read (ninp(0),'(a8)',err=91) id
      call xertst (streql(id,'*ttini'),
     &  'faulty input id instead of *ttini, id='//id)
      read (ninp(0),*,err=92) ttini
      call xertst (0.lt.ttini, 'faulty input ttini')
*   ..read and echo code internal parameters                            !srv 22.12.09 {
      write (nout(0),'(/2x,a)') 'non-default internal parameters:'
    1 continue
       read (ninp(0),'(a)',end=2,err=94) line
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

*   ..read the header of geometry file
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
*     (isymm is now stored in the b2cmfs common)               ! xpb 10.08.00
*   ..read the mapping and geometry
      call read_b2fgmtry (ninp(1), mpg, geo)

*   ..initialise time
      tim = 0.0_R8

      if(read_b2fstate_file.ne.0) then
*  .. copy the data of first n=min(ns0,ns1) species
*     from the plasma state file b2fstate for ns0 species
*     into the initial plasma state file b2fstati for ns1 species
       filename='b2fstate'
       call find_file(filename,ex)
       if (ex) then
         write(*,*) ' Opening '//trim(filename)
         call cfopen (ninp(2),filename,'old','un*formatted')
       else
         call xerrab ('Could not find old b2fstate file')
       end if
*   ..obtain version number
       call cfverr (ninp(2), b2fstate_version)
       if (b2fstate_version .lt. '03.002.000') then
         call xerrab('Using old input format b2fstate. Run b2us first.')
       endif
*   ..obtain dimensions from plasma state file
       call cfruin (ninp(2), 3, idum, 'nCv,nFc,ns')
       write(*,*) ' nx0, ny0, ns0', idum(0:2)
       ns0 = idum(2)
       call xertst (idum(0).eq.nCv.and.idum(1).eq.nFc,
     &   'faulty input nCv, nFc from plasma state file')
       call xertst (1.le.ns0,
     &  'faulty input ns from plasma state file')
*    ..read header for plasma initial state
*      (zamin, zamax, zn, am from ninp(2) will be checked against the data
*      already read in)
       allocate(zcdum(0:ns0-1,0:4))
       call cfruch (ninp(2), 120, lblps, 'label')
       call b2ruzd (ninp(2), b2fstate_version, ns0,
     &   zcdum(0,1), zcdum(0,2), zcdum(0,3), zcdum(0,4), .false.)
       do is = 0, min(ns,ns0)-1
        call xertst (ltst(zcdum(is,1),zamin(is)),
     &  'faulty data for zamin')
        call xertst (ltst(zcdum(is,2),zamax(is)),
     &  'faulty data for zamax')
        call xertst (ltst(zcdum(is,3),zn(is)), 'faulty data for zn')
        call xertst (ltst(zcdum(is,4),am(is)), 'faulty data for am')
       enddo
*    ..create the initial state
       call createB2State (nCv, nFc, ns, st1)
*    ..read initial state from the file
       call read_b2fstate (ninp(2), nCv, nFc, ns0, st1)
       close(ninp(2))
       write(*,*) ' Closed b2fstate'
*   ..test plasma state and fluxes
       call b2xvps (nCv, nFc, ns0, st1%pl, st1%dv)
       call b2xvsg (nCv*ns0, st1%dv%kinrgy, 1, 'kinrgy', '.ge.')
*   ..initialise na, ua, kinrgy, uadia, ne
       do is = 0, min(ns,ns0)-1
         st%pl%na(:,is)         = st1%pl%na(:,is)
         st%pl%ua(:,is)         = st1%pl%ua(:,is)
         st%dv%kinrgy(:,is)     = st1%dv%kinrgy(:,is)
         st%dv%uadia(:,:,is)  = st1%dv%uadia(:,:,is)
         st%dv%vadia(:,:,is)  = st1%dv%vadia(:,:,is)
         st%dv%wadia(:,:,is)  = st1%dv%wadia(:,:,is)
         st%dv%vaecrb(:,:,is) = st1%dv%vaecrb(:,:,is)
       enddo
       do is = min(ns,ns0), ns-1
        call sfill (nCv, naini(is), st%pl % na(1,is), 1)
        call sfill (nCv, 0.0_R8, st%pl%ua(1,is), 1)
        call sfill (nCv, 0.0_R8, st%dv%kinrgy(1,is), 1)
        call sfill (nFc*2, 0.0_R8, st%dv%uadia(1,0,is), 1)
        call sfill (nFc*2, 0.0_R8, st%dv%vadia(1,0,is), 1)
        call sfill (nFc*2, 0.0_R8, st%dv%wadia(1,0,is), 1)
        call sfill (nFc*2, 0.0_R8, st%dv%vaecrb(1,0,is), 1)
       enddo
       do is = 0, ns-1
        call sfill (nCv, (zamin(is)+zamax(is))/2.0_R8, 
     .    st%rt%rza(1,is), 1)
       enddo
       call b2xpne (nCv,ns,st%rt%rza,st%pl%na,zero,st%dv%ne)
*   ..initialise te, ti, tn, po, kt, zt
       st%pl%te = st1%pl%te
       st%pl%ti = st1%pl%ti
       st%pl%tn = st1%pl%tn
       st%pl%po = st1%pl%po
       st%pl%kt = st1%pl%kt
       st%pl%zt = st1%pl%zt
*   ..initialise fluxes
       do is = 0, min(ns,ns0)-1
         st%dv%fna(:,:,is)         = st1%dv%fna(:,:,is)
         st%dv%fna_mdf(:,:,is)     = st1%dv%fna_mdf(:,:,is)
         st%dv%fna_fcor(:,:,is)    = st1%dv%fna_fcor(:,:,is)
         st%dv%fna_nodrift(:,:,is) = st1%dv%fna_nodrift(:,:,is)
         st%dv%fna_he(:,:,is)      = st1%dv%fna_he(:,:,is)
         st%dv%fnaPSch(:,:,is)     = st1%dv%fnaPSch(:,:,is)
         st%dv%fna_eir(:,:,is)     = st1%dv%fna_eir(:,:,is)
         st%dv%fna_32(:,:,is)      = st1%dv%fna_32(:,:,is)
         st%dv%fna_52(:,:,is)      = st1%dv%fna_52(:,:,is)
       enddo
       do is = min(ns,ns0), ns-1
         call sfill (nFc*2, 0.0_R8, st%dv%fna(1,0,is), 1)
         call sfill (nFc*2, 0.0_R8, st%dv%fna_mdf(1,0,is), 1)
         call sfill (nFc*2, 0.0_R8, st%dv%fna_fcor(1,0,is), 1)
         call sfill (nFc*2, 0.0_R8, st%dv%fna_nodrift(1,0,is), 1)
         call sfill (nFc*2, 0.0_R8, st%dv%fna_he(1,0,is), 1)
         call sfill (nFc*2, 0.0_R8, st%dv%fnaPSch(1,0,is), 1)
         call sfill (nFc*2, 0.0_R8, st%dv%fna_eir(1,0,is), 1)
         call sfill (nFc*2, 0.0_R8, st%dv%fna_32(1,0,is), 1)
         call sfill (nFc*2, 0.0_R8, st%dv%fna_52(1,0,is), 1)
       enddo
       st%dv%vedia     = st1%dv%vedia
       st%dv%veecrb    = st1%dv%veecrb
       st%dv%fne_eir   = st1%dv%fne_eir
       st%dv%fne_32    = st1%dv%fne_32
       st%dv%fne_52    = st1%dv%fne_52
       st%dv%fni_32    = st1%dv%fni_32
       st%dv%fni_52    = st1%dv%fni_52
       st%dv%fhe       = st1%dv%fhe
       st%dv%fhe_mdf   = st1%dv%fhe_mdf
       st%dv%fhePSch   = st1%dv%fhePSch
       st%dv%fhe_eir   = st1%dv%fhe_eir
       st%dv%fhi       = st1%dv%fhi
       st%dv%fhn       = st1%dv%fhn
       st%dv%fhi_mdf   = st1%dv%fhi_mdf
       st%dv%fhiPSch   = st1%dv%fhiPSch
       st%dv%fhi_eir   = st1%dv%fhi_eir
       st%dv%fch       = st1%dv%fch
       st%dv%fch_32    = st1%dv%fch_32
       st%dv%fch_52    = st1%dv%fch_52
       st%dv%fch_p     = st1%dv%fch_p
       st%dv%fchdia    = st1%dv%fchdia
       st%dv%fchin     = st1%dv%fchin
       st%dv%fchvispar = st1%dv%fchvispar
       st%dv%fchvisper = st1%dv%fchvisper
       st%dv%fchvisq   = st1%dv%fchvisq
       st%dv%fchinert  = st1%dv%fchinert
       st%dv%floe_noc  = st1%dv%floe_noc
       st%dv%floi_noc  = st1%dv%floi_noc
       call destroyB2State (st1) 
      else                                                              !srv 22.12.09 }
*   ..initialise na, ua, kinrgy, uadia, ne
       do is = 0, ns-1
        call sfill (nCv, naini(is), st%pl%na(1,is), 1)
        call sfill (nCv, 0.0_R8, st%pl%ua(1,is), 1)
        call sfill (nCV, 0.0_R8, st%dv%kinrgy(1,is), 1)
        call sfill (nFc*2, 0.0_R8, st%dv%uadia(1,0,is), 1)
        call sfill (nFc*2, 0.0_R8, st%dv%vadia(1,0,is), 1)
        call sfill (nFc*2, 0.0_R8, st%dv%wadia(1,0,is), 1)
        call sfill (nFc*2, 0.0_R8, st%dv%vaecrb(1,0,is), 1)
        call sfill (nCv, (zamin(is)+zamax(is))/2.0_R8,
     .    st%rt%rza(1,is), 1)
       enddo
       call b2xpne (nCv,ns,st%rt%rza,st%pl%na,zero,st%dv%ne)
       call sfill (nFc*2, 0.0_R8, st%dv%veecrb, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%vedia, 1)
*   ..initialise te, ti, tn, po, kt, zt
       call sfill (nCv, ttini*ev, st%pl%te, 1)
       call sfill (nCv, ttini*ev, st%pl%ti, 1)
       call sfill (nCv, ttini*ev, st%pl%tn, 1)
       call sfill (nCv, 3.1_R8*ttini, st%pl%po, 1)
       call sfill (nCv, 0.0_R8, st%pl%kt, 1)
       call sfill (nCv, 0.0_R8, st%pl%zt, 1)
*   ..initialise fluxes
       call sfill (nFc*2*ns, 0.0_R8, st%dv%fna, 1)
       call sfill (nFc*2*ns, 0.0_R8, st%dv%fna_mdf, 1)
       call sfill (nFc*2*ns, 0.0_R8, st%dv%fna_fcor, 1)
       call sfill (nFc*2*ns, 0.0_R8, st%dv%fna_nodrift, 1)
       call sfill (nFc*2*ns, 0.0_R8, st%dv%fna_he, 1)
       call sfill (nFc*2*ns, 0.0_R8, st%dv%fnaPSch, 1)
       call sfill (nFc*2*ns, 0.0_R8, st%dv%fna_eir, 1)
       call sfill (nFc*2*ns, 0.0_R8, st%dv%fna_32, 1)
       call sfill (nFc*2*ns, 0.0_R8, st%dv%fna_52, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fne_eir, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fhe_mdf, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fhe_eir, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fhePSch, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fhi_mdf, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fhi_eir, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fhiPSch, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fhe, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fhi, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fhn, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fch, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fne_32, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fne_52, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fni_32, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fni_52, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fch_32, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fch_52, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fch_p, 1)                       !srv 10.06.08
       call sfill (nFc*2, 0.0_R8, st%dv%fchdia, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fchin, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fchvispar, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fchvisper, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fchvisq, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%fchinert, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%floe_noc, 1)
       call sfill (nFc*2, 0.0_R8, st%dv%floi_noc, 1)
      endif                                                             !srv 22.12.09

*   ..write b2fstati
      idum(0) = nCv
      idum(1) = nFc
      idum(2) = ns
      call cfverw (nout(1), newversion)
      call cfwuin (nout(1), 3, idum, 'nCv,nFc,ns')
      call cfwuch (nout(1), 120, lblps, 'label')
      call b2wuzd (nout(1), newversion, ns, zamin, zamax, zn, am)
      call write_b2fstate (nout(1), nCv, nFc, ns, st)

*   ..for underlying structured grids, convert to structured output and write b2fstate_st
      if (mpg%isClassicalGrid.eq.1) then
       call to_struct_plasma (ns,mpg,st)
      endif

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
