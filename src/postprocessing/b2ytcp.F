*-----------------------------------------------------------------------
*.specification

      subroutine b2ytcp (nget, nput, ns, ns1,
     * zamin1, zamax1, zn1, am1, is_neutral1, bundle_frac)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      implicit none
*   ..input arguments (unchanged on exit)
      integer nget, nput, ns, ns1
      real (kind=R8) ::
     *  zamin1(0:ns1-1), zamax1(0:ns1-1), zn1(0:ns1-1), am1(0:ns1-1),
     *  bundle_frac(0:ns1-1,0:ns-1)
      logical is_neutral1(0:ns1-1)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YTCP reads from a formatted file some parameters in the
*     physics common blocks and rewrites them after conversion.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine reads from unit (nget) the values of some parameters
*     in the physics common blocks: /b2cmpa/, /b2cmpb/, and /b2cmpt/.
*     Unit nget must be connected to a formatted file.
*     The input file must contain a sequence of segments (the sequence
*     may be empty). Each segment begins with a record that starts with
*     a text of the form '*heading' or '*HEADING', in which heading
*     identifies the set of parameters for which values are provided in
*     the subsequent records of that segment. The subsequent records
*     are generally in free format; see the code to discern the precise
*     requirements.
*
*     Assignments are made only to those elements of the common blocks
*     for which values are provided in the input file. The other
*     elements are left unchanged.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nget - integer, input.
*     nget specifies a unit number for formatted input.
*
*     nput - integer, input.
*     nput specifies a unit number for formatted output.
*
*     ns - integer, input.
*     ns specifies the number of atomic species in the input set.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*     ns1 - integer, input.
*     ns1 specifies the number of atomic species in the output set.
*     The species are indexed by (0:ns1-1).
*     It will hold that 1.le.ns1.
*
*     (parameters in common:)
*
*     /b2cmpa/ - input/output.
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
*     /b2cmpb/ - input/output.
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
*     /b2cmpt/ - input/output.
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
*         cfsig(3)*te(,)/ev/abs(bb(,,3))/16+
*         cfsig(4)*exp(-cfsig(5)*(y/L)**2)
*     that specifies the diffusivity corresponding to the anomalous
*     electrical conductivity, which is
*       sig(,) = qe*ne(,)*dfsig(,).
*     L is the local radial width of the plasma and
*     y is the radial distance from the inner (South) boundary.
*     If cfsig(7) is not zero, the sig0 array is flux-scaled with
*     a scaling factor of cfsig(7)*(hy1(,)/hy1(ixref,))**2
*     where hy1 is the radial width of the cell and ixref is a reference
*     poloidal position (usually the outer midplane or 'jxa').
*     cfsig(2) and (6) are not used in the default version of the code.
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
*     cflim(3)=cflal specifies the flux limit factor for the
*     thermo-electric coefficient.
*     Zero settings turn off these flux limits.
*     cflim(4:7) are not used in the default version of the code.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, is0, i, ir0, ireg, mode, nsm
      real (kind=R8) :: zpdum(0:ns-1,0:3)
      real (kind=R8) :: cftmp(0:nsdecl-1)
      character id*8, id0*8, id1*72
      character*(*) fmt0, fmt1
      parameter (fmt0='(0p,1x,'''''''',i3,'''''''',1x,4(f9.3,1x))')
      parameter (fmt1='(1p,1x,'''''''',i3,'''''''',1x,8(e9.2,1x))')
*     (The sequence of eight quotes is not a joke. Because it occurs
*     inside a character string it represents four quotes, which in
*     turn represents one quoted quote in the output format.)
      logical oldversion
*   ..procedures
      logical streql
      external subini, subend, xertst, xerrab, streql
      intrinsic max

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     mode - integer.
*     mode indicates what group of data is being read: mode.eq.0 for
*     /b2cmpa/, mode.eq.1 for /b2cmpb/, and mode.eq.2 for /b2cmpt/.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ytcp')
*   ..test nget, nput
      call xertst (1.le.nget, 'faulty argument nget')
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test ns, ns1
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (1.le.ns1, 'faulty argument ns1')
*   ..test nsdecl
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

* ..read parameters
*   ..initialise mode
      mode = -1
*   ..initialise nsm
      nsm = max(ns,ns1)
*   ..start of loop
    1 continue
*   ..read id
      read (nget,'(a8)',err=91,end=2) id
*   ..test and branch
      if (streql(id,'*endphy')) then
       write (nput,'(a)') '*endphy'
*    ..end of sequence
       goto 2
      else if (streql(id,'*b2cmpa')) then
       mode = 0
*   ..write heading for /b2cmpa/
       write (nput,'(a)')
     &  '*b2cmpa         basic parameters'
      else if (streql(id,'*b2cmpb')) then
       mode = 1
       ireg = -1
*   ..write heading for /b2cmpb/
       write (nput,'(a)')
     &  '*b2cmpb         boundary conditions'
      else if (streql(id,'*b2cmpt')) then
       mode = 2
*   ..write heading for /b2cmpt/
       write (nput,'(a)')
     &  '*b2cmpt         transport coefficients'
      else if (streql(id,'*specs')) then
*    ..test mode
       call xertst (mode.eq.0, 'unexpected header '//id)
*    ..read species data
*    ..Two input formats are accepted:
*    ..old version: za, zn, am, za2
*    ..new version: zamin, zamax, zn, am
*    ..The output format is new version only
       oldversion = .true.
       do is = 0, ns-1
        read (nget,*,err=92) id0, (zpdum(is,i),i=0,3)
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
*   ..write species data
       write (nput,'(a)')
     &  '*specs      zamin    zamax    zn      am'
       write (nput,fmt0) 
     &  (is,zamin1(is),zamax1(is),zn1(is),am1(is),is=0,ns1-1)
      else if (streql(id,'*cbregs')) then
*    ..test mode
       call xertst (mode.eq.1, 'unexpected header '//id)
*    ..read region specs
       cbirso = 0
       read (nget,*,err=92) id0, cbnrso,
     &   (cbrbrk(ir0),ir0=cbirso,cbirso+cbnrso-1)
       call xertst (streql(id0,'south').and.1.le.cbnrso.and.
     &   cbirso+cbnrso.le.cbnrmx, 'faulty input under *cbregs')
       cbirno = cbirso+cbnrso
       read (nget,*,err=92) id0, cbnrno,
     &   (cbrbrk(ir0),ir0=cbirno,cbirno+cbnrno-1)
       call xertst (streql(id0,'north').and.1.le.cbnrno.and.
     &   cbirno+cbnrno.le.cbnrmx, 'faulty input under *cbregs')
       cbirwe = cbirno+cbnrno
       read (nget,*,err=92) id0, cbnrwe,
     &   (cbrbrk(ir0),ir0=cbirwe,cbirwe+cbnrwe-1)
       call xertst (streql(id0,'west').and.1.le.cbnrwe.and.
     &   cbirwe+cbnrwe.le.cbnrmx, 'faulty input under *cbregs')
       cbirea = cbirwe+cbnrwe
       read (nget,*,err=92) id0, cbnrea,
     &   (cbrbrk(ir0),ir0=cbirea,cbirea+cbnrea-1)
       call xertst (streql(id0,'east').and.1.le.cbnrea.and.
     &   cbirea+cbnrea.le.cbnrmx, 'faulty input under *cbregs')
*   ..write region specs
       write (nput,'(a)')
     &  '*cbregs   nr    breakpoints'
       write (nput,'(1x,a,t12,i2,(t17,7f9.3))')
     &  '''south''', cbnrso, (cbrbrk(ireg),ireg=cbirso,cbirso+cbnrso-1)
       write (nput,'(1x,a,t12,i2,(t17,7f9.3))')
     &  '''north''', cbnrno, (cbrbrk(ireg),ireg=cbirno,cbirno+cbnrno-1)
       write (nput,'(1x,a,t12,i2,(t17,7f9.3))')
     &  '''west''', cbnrwe, (cbrbrk(ireg),ireg=cbirwe,cbirwe+cbnrwe-1)
       write (nput,'(1x,a,t12,i2,(t17,7f9.3))')
     &  '''east''', cbnrea, (cbrbrk(ireg),ireg=cbirea,cbirea+cbnrea-1)
      else if (streql(id,'*region')) then
*    ..test mode
       call xertst (mode.eq.1, 'unexpected header '//id)
*    ..read ireg
       read (nget,*,err=92) ireg, id1
       if (cbirso.le.ireg.and.ireg.lt.cbirso+cbnrso) then
        call xertst (streql(id1(1:5),'south'),
     &    'faulty input under *region: '//trim(id1))
        write (nput,'(a/2x,i2,2x,a)') '*region', ireg, 
     &    ''''//trim(id1)//''''
       else if (cbirno.le.ireg.and.ireg.lt.cbirno+cbnrno) then
        call xertst (streql(id1(1:5),'north'),
     &    'faulty input under *region: '//trim(id1))
        write (nput,'(a/2x,i2,2x,a)') '*region', ireg, 
     &    ''''//trim(id1)//''''
       else if (cbirwe.le.ireg.and.ireg.lt.cbirwe+cbnrwe) then
        call xertst (streql(id1(1:4),'west'),
     &    'faulty input under *region: '//trim(id1))
        write (nput,'(a/2x,i2,2x,a)') '*region', ireg,
     &    ''''//trim(id1)//''''
       else if (cbirea.le.ireg.and.ireg.lt.cbirea+cbnrea) then
        call xertst (streql(id1(1:4),'east'),
     &    'faulty input under *region: '//trim(id1))
        write (nput,'(a/2x,i2,2x,a)') '*region', ireg,
     &    ''''//trim(id1)//''''
       else
        call xerrab ('faulty input under *region: out-of-range')
       endif
      else if (streql(id,'*cbsna')) then
*    ..test mode and ireg
       call xertst (mode.eq.1, 'unexpected header '//id)
       call xertst (0.le.ireg, 'unexpected header '//id//
     &   ', region is not defined')
*    ..read cbsna
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cbsna(i,is,ireg),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cbsna')
       enddo
*    ..convert cbsna
       do i = 0, 7
         cftmp(0:ns-1) = cbsna(i,0:ns-1,ireg)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cbsna(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
       enddo
*    ..write cbsna
       write (nput,'(a)')
     &   '*cbsna     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cbsna(i,is,ireg),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cbsmo')) then
*    ..test mode and ireg
       call xertst (mode.eq.1, 'unexpected header '//id)
       call xertst (0.le.ireg, 'unexpected header '//id//
     &   ', region is not defined')
*    ..read cbsmo
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cbsmo(i,is,ireg),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cbsmo')
       enddo
*    ..convert cbsmo
       do i = 0, 7
         cftmp(0:ns-1) = cbsmo(i,0:ns-1,ireg)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cbsmo(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
       enddo
*    ..write cbsmo
       write (nput,'(a)')
     &   '*cbsmo     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cbsmo(i,is,ireg),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cbshi')) then
*    ..test mode and ireg
       call xertst (mode.eq.1, 'unexpected header '//id)
       call xertst (0.le.ireg, 'unexpected header '//id//
     &   ', region is not defined')
*    ..read cbshi
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cbshi(i,is,ireg),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cbshi')
       enddo
*    ..convert cbshi
       do i = 0, 7
         cftmp(0:ns-1) = cbshi(i,0:ns-1,ireg)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cbshi(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
       enddo
*    ..write cbshi
       write (nput,'(a)')
     &   '*cbshi     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cbshi(i,is,ireg),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cbshe')) then
*    ..test mode and ireg
       call xertst (mode.eq.1, 'unexpected header '//id)
       call xertst (0.le.ireg, 'unexpected header '//id//
     &   ', region is not defined')
*    ..read cbshe
       read (nget,*,err=92) id0, (cbshe(i,ireg),i=0,7)
       read (id0,'(i8)',err=93) is0
       call xertst (is0.eq.-1, 'faulty input under *cbshe')
*    ..write cbshe
       write (nput,'(a)')
     &   '*cbshe     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) -1, (cbshe(i,ireg),i=0,7)
      else if (streql(id,'*cbsch')) then
*    ..test mode and ireg
       call xertst (mode.eq.1, 'unexpected header '//id)
       call xertst (0.le.ireg, 'unexpected header '//id//
     &   ', region is not defined')
*    ..read cbsch
       read (nget,*,err=92) id0, (cbsch(i,ireg),i=0,7)
       read (id0,'(i8)',err=93) is0
       call xertst (is0.eq.-1, 'faulty input under *cbsch')
*    ..write cbsch
       write (nput,'(a)')
     &   '*cbsch     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) -1, (cbsch(i,ireg),i=0,7)
      else if (streql(id,'*cbrec')) then
*    ..test mode and ireg
       call xertst (mode.eq.1, 'unexpected header '//id)
       call xertst (0.le.ireg, 'unexpected header '//id//
     &   ', region is not defined')
*    ..read cbrec
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cbrec(i,is,ireg),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cbrec')
       enddo
*    ..convert cbrec
       do i = 0, 7
         cftmp(0:ns-1) = cbrec(i,0:ns-1,ireg)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cbrec(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
       enddo
*    ..write cbrec
       write (nput,'(a)')
     &   '*cbrec     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cbrec(i,is,ireg),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cbmsa')) then
*    ..test mode and ireg
       call xertst (mode.eq.1, 'unexpected header '//id)
       call xertst (0.le.ireg, 'unexpected header '//id//
     &   ', region is not defined')
*    ..read cbmsa
       read (nget,*,err=92) id0, (cbmsa(i,ireg),i=0,7)
       read (id0,'(i8)',err=93) is0
       call xertst (is0.eq.-1, 'faulty input under *cbmsa')
*    ..write cbmsa
       write (nput,'(a)')
     &   '*cbmsa     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) -1, (cbmsa(i,ireg),i=0,7)
      else if (streql(id,'*cbmsb')) then
*    ..test mode and ireg
       call xertst (mode.eq.1, 'unexpected header '//id)
       call xertst (0.le.ireg, 'unexpected header '//id//
     &   ', region is not defined')
*    ..read cbmsb
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cbmsb(i,is,ireg),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cbmsb')
       enddo
*    ..convert cbmsb
       do i = 0, 7
         cftmp(0:ns-1) = cbmsb(i,0:ns-1,ireg)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cbmsb(i,0:ns1-1,ireg) = cftmp(0:ns1-1)
       enddo
*    ..write cbmsb
       write (nput,'(a)')
     &   '*cbmsb     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cbmsb(i,is,ireg),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cfdf0')) then
*    ..test mode
       call xertst (mode.eq.2, 'unexpected header '//id)
*    ..read cfdf0
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cfdf0(i,is),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cfdf0')
       enddo
*    ..convert cfdf0
       do i = 0, 7
         cftmp(0:ns-1) = cfdf0(i,0:ns-1)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cfdf0(i,0:ns1-1) = cftmp(0:ns1-1)
       enddo
*    ..write cfdf0
       write (nput,'(a)')
     &   '*cfdf0     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cfdf0(i,is),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cfdna')) then
*    ..test mode
       call xertst (mode.eq.2, 'unexpected header '//id)
*    ..read cfdna
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cfdna(i,is),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cfdna')
       enddo
*    ..convert cfdna
       do i = 0, 7
         cftmp(0:ns-1) = cfdna(i,0:ns-1)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cfdna(i,0:ns1-1) = cftmp(0:ns1-1)
       enddo
*    ..write cfdna
       write (nput,'(a)')
     &   '*cfdna     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cfdna(i,is),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cfdpa')) then
*    ..test mode
       call xertst (mode.eq.2, 'unexpected header '//id)
*    ..read cfdpa
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cfdpa(i,is),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cfdpa')
       enddo
*    ..convert cfdpa
       do i = 0, 7
         cftmp(0:ns-1) = cfdpa(i,0:ns-1)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cfdpa(i,0:ns1-1) = cftmp(0:ns1-1)
       enddo
*    ..write cfdpa
       write (nput,'(a)')
     &   '*cfdpa     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cfdpa(i,is),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cfvla')) then
*    ..test mode
       call xertst (mode.eq.2, 'unexpected header '//id)
*    ..read cfvla
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cfvla(i,is),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cfvla')
       enddo
*    ..convert cfvla
       do i = 0, 7
         cftmp(0:ns-1) = cfvla(i,0:ns-1)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cfvla(i,0:ns1-1) = cftmp(0:ns1-1)
       enddo
*    ..write cfvla
       write (nput,'(a)')
     &   '*cfvla     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cfvla(i,is),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cfvsa')) then
*    ..test mode
       call xertst (mode.eq.2, 'unexpected header '//id)
*    ..read cfvsa
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cfvsa(i,is),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cfvsa')
       enddo
*    ..convert cfvsa
       do i = 0, 7
         cftmp(0:ns-1) = cfvsa(i,0:ns-1)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cfvsa(i,0:ns1-1) = cftmp(0:ns1-1)
       enddo
*    ..write cfvsa
       write (nput,'(a)')
     &   '*cfvsa     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cfvsa(i,is),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cfhci')) then
*    ..test mode
       call xertst (mode.eq.2, 'unexpected header '//id)
*    ..read cfhci
       do is = 0, ns-1
        read (nget,*,err=92) id0, (cfhci(i,is),i=0,7)
        read (id0,'(i8)',err=93) is0
        call xertst (is0.eq.is, 'faulty input under *cfhci')
       enddo
*    ..convert cfhci
       do i = 0, 7
         cftmp(0:ns-1) = cfhci(i,0:ns-1)
         call convert_bundle_r
     .    (ns1,ns,nsm,is_neutral1,bundle_frac,cftmp)
         cfhci(i,0:ns1-1) = cftmp(0:ns1-1)
       enddo
*    ..write cfhci
       write (nput,'(a)')
     &   '*cfhci     (0)       (1)       (2)       (3)       (4)'//
     &   '       (5)       (6)       (7)'
       write (nput,fmt1) (is,(cfhci(i,is),i=0,7),is=0,ns1-1)
      else if (streql(id,'*cfhce')) then
*    ..test mode
       call xertst (mode.eq.2, 'unexpected header '//id)
*    ..read cfhce
       read (nget,*,err=92) id0, (cfhce(i),i=0,7)
       read (id0,'(i8)',err=93) is0
       call xertst (is0.eq.-1, 'faulty input under *cfhce')
*    ..write cfhce
       write (nput,'(a)')
     &  '*cfhce     (0)       (1)       (2)       (3)       (4)'//
     &  '       (5)       (6)       (7)'
       write (nput,fmt1) -1, (cfhce(i),i=0,7)
      else if (streql(id,'*cfsig')) then
*    ..test mode
       call xertst (mode.eq.2, 'unexpected header '//id)
*    ..read cfsig
       read (nget,*,err=92) id0, (cfsig(i),i=0,7)
       read (id0,'(i8)',err=93) is0
       call xertst (is0.eq.-1, 'faulty input under *cfsig')
*    ..write cfsig
       write (nput,'(a)')
     &  '*cfsig     (0)       (1)       (2)       (3)       (4)'//
     &  '       (5)       (6)       (7)'
       write (nput,fmt1) -1, (cfsig(i),i=0,7)
      else if (streql(id,'*cfalf')) then
*    ..test mode
       call xertst (mode.eq.2, 'unexpected header '//id)
*    ..read cfalf
       read (nget,*,err=92) id0, (cfalf(i),i=0,7)
       read (id0,'(i8)',err=93) is0
       call xertst (is0.eq.-1, 'faulty input under *cfalf')
*    ..write cfalf
       write (nput,'(a)')
     &  '*cfalf     (0)       (1)       (2)       (3)       (4)'//
     &  '       (5)       (6)       (7)'
       write (nput,fmt1) -1, (cfalf(i),i=0,7)
      else if (streql(id,'*cflim')) then
*    ..test mode
       call xertst (mode.eq.2, 'unexpected header '//id)
*    ..read cflim
       read (nget,*,err=92) id0, (cflim(i),i=0,7)
       read (id0,'(i8)',err=93) is0
       call xertst (is0.eq.-1, 'faulty input under *cflim')
*    ..write cflim
       write (nput,'(a)')
     &  '*cflim     (0)       (1)       (2)       (3)       (4)'//
     &  '       (5)       (6)       (7)'
       write (nput,fmt1) -1, (cflim(i),i=0,7)
      else
       call xerrab ('faulty input of id, id='//id)
      endif
*   ..loop
      goto 1
*   ..end-of-file exit
    2 continue

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read header id')
   92 call xerrab ('error trying to read after '//id)
   93 call xerrab ('error trying to read species index after '//id)

*-----------------------------------------------------------------------
*.end b2ytcp

      end subroutine b2ytcp
