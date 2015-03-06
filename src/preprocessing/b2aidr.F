*-----------------------------------------------------------------------
*.specification

      subroutine b2aidr (ninp, nout, nx, ny, ns)
      use b2mod_types
      use b2mod_rates
      use b2mod_geo
      use b2mod_time
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
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
*     state for b2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine prepares and writes to an external file an initial
*     plasma state for b2.
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
      integer n2, is, i, is0, idum(0:9), zero(-1:nx,-1:ny)
      integer, save :: read_b2fstate = 0                              !srv 22.12.09
      real (kind=R8) :: t0, t1                                        !srv 22.12.09
      real (kind=R8), allocatable :: zcdum(:,:), 
     &  na1(:,:,:), ua1(:,:,:),  ne1(:,:), uadia1(:,:,:,:), 
     &  te1(:,:), ti1(:,:), po1(:,:), fna1(:,:,:,:), fhe1(:,:,:),
     &  fhi1(:,:,:), fch1(:,:,:), fch_321(:,:,:), fch_521(:,:,:),
     &  fch_p1(:,:,:), kinrgy1(:,:,:)   
      real (kind=R8) :: zpdum(0:ns-1,0:3)
      character lblps*120, lblgm*120, id*8, id0*8
      character idout*80, cnamip*80, cvalip*80                        !srv 22.12.09
      real (kind=R8) ::
     *  naini(0:ns-1), ttini
      logical oldversion
*   ..procedures
      logical streql
      logical ltst                                                    !srv 22.12.09
      external subini, subend, xertst, xerrab, sfill, streql,
     &  cfwuch, cfwuin
      external b2xvsg, b2xxid, b2rflb, b2wups, b2rugm, b2xpne, b2wuzd,
     &  b2ruzd, b2rups
      ltst(t0,t1) = abs(t0-t1).le.1.0e-6_R8*(abs(t0)+abs(t1))         !srv 22.12.09

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
*
*     uadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), uadia(ix,iy,0,is) specifies
*     the diamagnetic drift velocity of species (is) on the (ix,iy) cell, and
*     uadia(ix,iy,1,is) the radial drift velocity.
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
*     heat flux through the face between the (ix,iy-1) and (ix,iy)
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
     &  1.le.nout(0).and.1.le.nout(1),
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
*     (ttini is provided in units of ev)
      read (ninp(1),'(a8)',err=91) id
      call xertst (streql(id,'*ttini'),
     &  'faulty input id instead of *ttini, id='//id)
      read (ninp(1),*,err=92) ttini
      call xertst (0.lt.ttini, 'faulty input ttini')
*   ..read and echo code internal parameters                      !srv 22.12.09 {
      write (nout(0),'(/2x,a)') 'non-default internal parameters:'
    1 continue
       read (ninp(1),*,end=2,err=94) cnamip, cvalip
       write (nout(0),'(4x,a,2x,a)') cnamip, cvalip
       call ipsetc (cnamip, cvalip)
       goto 1
    2 continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'
      call ipgeti ('b2aidr_read_b2fstate', read_b2fstate)

      if(read_b2fstate.ne.0) then
*  .. copy the data of first n=min(ns0,ns1) species 
*     from the plasma state file b2fstate for ns0 species
*     into the initial plasma state file b2fstati for ns1 species 
       write(*,*) ' Opening b2fstate'
       call cfopen (ninp(3),'b2fstate','old','un*formatted')
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
     &  zcdum(0,1), zcdum(0,2), zcdum(0,3), zcdum(0,4))
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
     &  ne1(-1:nx,-1:ny), uadia1(-1:nx,-1:ny,0:1,0:idum(2)-1), 
     &  te1(-1:nx,-1:ny), ti1(-1:nx,-1:ny), po1(-1:nx,-1:ny),
     &  fna1(-1:nx,-1:ny,0:1,0:idum(2)-1), fhe1(-1:nx,-1:ny,0:1),
     &  fhi1(-1:nx,-1:ny,0:1), fch1(-1:nx,-1:ny,0:1),
     1  fch_321(-1:nx,-1:ny,0:1), fch_521(-1:nx,-1:ny,0:1),
     2  fch_p1(-1:nx,-1:ny,0:1), kinrgy1(-1:nx,-1:ny,0:idum(2)-1) )
*   ..read initial state from the file
       call b2rups (ninp(3), idum(0), idum(1), idum(2), ne1, na1, ua1, 
     &  uadia1, te1, ti1, po1, fna1, fhe1, fhi1, fch1, fch_321, fch_521,  
     &  fch_p1, kinrgy1, b2fstate_version)
       close(ninp(3))
       write(*,*) ' Closed b2fstate'
*   ..test plasma state and fluxes
       call b2xvps (idum(0), idum(1), idum(2), 
     &  ne1, na1, ua1, te1, ti1, po1, fna1, fhe1, fhi1, fch1)
       call b2xvsg (n2*idum(2), kinrgy1, 1, 'kinrgy', '.ge.')
*   ..initialise na, ua, kinrgy, uadia, ne
       do is = 0, min(ns,idum(2))-1
	  na(:,:,is)     = na1(:,:,is)
	  ua(:,:,is)     = ua1(:,:,is)
	  kinrgy(:,:,is) = kinrgy1(:,:,is)
	  uadia(:,:,:,is)  = uadia1(:,:,:,is)
       enddo
       do is = min(ns,idum(2)), ns-1
        call sfill (n2, naini(is), na(-1,-1,is), 1)
        call sfill (n2, 0.0_R8, ua(-1,-1,is), 1)
        call sfill (n2, 0.0_R8, kinrgy(-1,-1,is), 1)
        call sfill (n2*2, 0.0_R8, uadia(-1,-1,0,is), 1)
       enddo
       do is = 0, ns-1
        call sfill (n2, (zamin(is)+zamax(is))/2.0_R8, rza(-1,-1,is), 1)
       enddo
       call b2xpne (nx, ny, ns, rza, na, zero, ne)
*   ..initialise te, ti, po
       te = te1
       ti = ti1
       po = po1
*   ..initialise fluxes
       do is = 0, min(ns,idum(2))-1
	  fna(:,:,:,is)     = fna1(:,:,:,is)
       enddo
       do is = min(ns,idum(2)), ns-1
        call sfill (n2*2, 0.0_R8, fna(-1,-1,0,is), 1)
       enddo
	 fhe    = fhe1
	 fhi    = fhi1
	 fch    = fch1
	 fch_32 = fch_321
	 fch_52 = fch_521
	 fch_p  = fch_p1
         deallocate( zcdum, na1, ua1, ne1, uadia1, te1, ti1, po1,
     &    fna1, fhe1, fhi1, fch1, fch_321, fch_521, fch_p1, kinrgy1)
      else                                                           !srv 22.12.09 }
*   ..initialise na, ua, kinrgy, uadia, ne
       do is = 0, ns-1
        call sfill (n2, naini(is), na(-1,-1,is), 1)
        call sfill (n2, 0.0_R8, ua(-1,-1,is), 1)
        call sfill (n2, 0.0_R8, kinrgy(-1,-1,is), 1)
        call sfill (n2*2, 0.0_R8, uadia(-1,-1,0,is), 1)
        call sfill (n2, (zamin(is)+zamax(is))/2.0_R8, rza(-1,-1,is), 1)
       enddo
       call b2xpne (nx, ny, ns, rza, na, zero, ne)
*   ..initialise te, ti, po
       call sfill (n2, ttini*ev, te, 1)
       call sfill (n2, ttini*ev, ti, 1)
       call sfill (n2, 3.1_R8*ttini, po, 1)
*   ..initialise fluxes
       call sfill (n2*2*ns, 0.0_R8, fna, 1)
       call sfill (n2*2, 0.0_R8, fhe, 1)
       call sfill (n2*2, 0.0_R8, fhi, 1)
       call sfill (n2*2, 0.0_R8, fch, 1)
       call sfill (n2*2, 0.0_R8, fch_32, 1)
       call sfill (n2*2, 0.0_R8, fch_52, 1)
       call sfill (n2*2, 0.0_R8, fch_p, 1)                           !srv 10.06.08
      endif                                                          !srv 22.12.09
*   ..initialise time
      tim = 0.0_R8
*   ..obtain version number						 !xpb
      call cfverr (ninp(2), b2fgmtry_version)				 !xpb
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(2), 2, idum, 'nx,ny')
      call xertst (idum(0).eq.nx.and.idum(1).eq.ny,
     &  'faulty input nx, ny from geometry file')
*   ..read header for geometry
      call cfruch (ninp(2), 120, lblgm, 'label')
      call cfruin (ninp(2), 1, idum, 'isymm')
      isymm = idum(0)
*     (isymm is now stored in the b2cmfs common)               ! xpb 10.08.00
*   ..read geometry
      call b2rugm (ninp(2), nx, ny, crx, cry, fpsi,
     &  ffbz, bb, vol, hx, hy, qz, qc, gs, pbs, wbbl, wbbr, wbbv)

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
      write (nout(0),'(6x,a3,2x,1pe11.2)') 'te:', ttini
      write (nout(0),'(6x,a3,2x,1pe11.2)') 'ti:', ttini
      write (nout(0),'(6x,a3,2x,1pe11.2)') 'po:', 3.1_R8*ttini
*   ..test plasma state
      call b2xvps (nx, ny, ns, 
     .  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
      call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
*   ..save dimensions
      idum(0) = nx
      idum(1) = ny
      idum(2) = ns
      call cfverw (nout(1), newversion)					 !xpb
      call cfwuin (nout(1), 3, idum, 'nx,ny,ns')
*   ..save label
      call cfwuch (nout(1), 120, lblps, 'label')
*   ..save /b2cmpa/
      call b2wuzd (nout(1), newversion, ns, zamin, zamax, zn, am)
*   ..save state
      call b2wups (nout(1), nx, ny, ns, ne, na, ua, uadia, te,		 !xpb
     &  ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy)

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
     &  ('  species',i3,',  atomic charge range from',f6.2,
     &   ' to',f6.2,
     &   ',  nuclear charge',f6.2,',  atomic mass',f6.2))

*-----------------------------------------------------------------------
*.end b2aidr

      end subroutine b2aidr
