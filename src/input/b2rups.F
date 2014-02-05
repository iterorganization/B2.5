*-----------------------------------------------------------------------
*.specification

      subroutine b2rups (nget, nx, ny, ns, 
     &  ne, na, ua, uadia, te, ti, po,
     &  fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy, 
     &  b2fstat_version)
      use b2mod_types
      use b2mod_rates
      use b2mod_time
      use b2mod_constants
      use b2mod_b2cmpa
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
     2  fch_p(-1:nx,-1:ny,0:1), kinrgy(-1:nx,-1:ny,0:ns-1)    !srv 10.06.08

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
*    
*     uadia - (-1:nx,-1:ny,0:1,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), uadia(ix,iy,0,is) specifies
*     the diamagnetic drift velocity of species (is) on the (ix,iy) cell, and
*     uadia(ix,iy,1,is) the radial drift velocity.
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
*     (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fnax(ix,-1:ny,0:ns-1) will hold 0.
*     For (ix,iy,is) in (-1:nx,0:ny,0:ns-1), fnay(ix,iy,is) specifies
*     the flux of atoms of species (is) through the face between the
*     (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fnay(-1:nx,iy,0:ns-1) will hold 0.
*
*     fhe - (-1:nx,-1:ny,0:1) real*8 array, output.
*     (Let fhex(,)=fhe(,,0) and fhey(,)=fhe(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhex(ix,iy) specifies the electron
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhex(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhey(ix,iy) specifies the electron
*     heat flux through the face between the (ix,iy-1) and (ix,iy)
*     cells. For iy.eq.-1, fhey(-1:nx,iy) will hold 0.
*
*     fhi - (-1:nx,-1:ny,0:1) real*8 array, output.
*     (Let fhix(,)=fhi(,,0) and fhiy(,)=fhi(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fhix(ix,iy) specifies the all atom
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For ix.eq.-1, fhix(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fhiy(ix,iy) specifies the all atom
*     heat flux through the face between the (ix-1,iy) and (ix,iy)
*     cells. For iy.eq.-1, fhiy(-1:nx,iy) will hold 0.
*
*     fch - (-1:nx,-1:ny,0:1) real*8 array, output.
*     (Let fchx(,)=fch(,,0) and fchy(,)=fch(,,1).)
*     For (ix,iy) in (0:nx,-1:ny), fchx(ix,iy) specifies the electric
*     current through the face between the (ix-1,iy) and (ix,iy) cells.
*     For ix.eq.-1, fchx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), fchy(ix,iy) specifies the electric
*     current through the face between the (ix,iy-1) and (ix,iy) cells.
*     For iy.eq.-1, fchy(-1:nx,iy) will hold 0.
*
*     kinrgy - (-1:nx,-1:ny,0:ns-1) real*8 array.
*     For (ix,iy,is) in (-1:nx,-1:ny,0:ns-1), kinrgy(ix,iy,is) specifies
*     the parallel kinetic energy of a particle of species (is) in 
*     cell (ix,iy).
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
      real (kind=R8) :: time(1)
*   ..procedures
      external subini, subend, xertst, cfrure
      external b2xvsg, b2xpne

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2rups')
*   ..set n2
      n2 = (nx+2)*(ny+2)
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
        call b2xpne (nx, ny, ns, rza, na, ne)
      endif
      call cfrure (nget, n2*ns, ua, 'ua')
      if(b2fstat_version.ge.'01.001.000') then                           !xpb
        call cfrure (nget, n2*2*ns, uadia, 'uadia')                      !xpb
      else                                                               !dpc
        uadia=0.0_R8                                                     !dpc
      endif                                                              !xpb
      call cfrure (nget, n2, te, 'te')
      call cfrure (nget, n2, ti, 'ti')
      call cfrure (nget, n2, po, 'po')
*   ..read fluxes
      call cfrure (nget, n2*2*ns, fna, 'fna')
      call cfrure (nget, n2*2, fhe, 'fhe')
      call cfrure (nget, n2*2, fhi, 'fhi')
      call cfrure (nget, n2*2, fch, 'fch')
      if(b2fstat_version.ge.'01.001.005') then                          !xpb
        call cfrure (nget, n2*2, fch_32, 'fch_32')
        call cfrure (nget, n2*2, fch_52, 'fch_52')
      else                                                              !dpc
        fch_32(:,:,0)=fch(:,:,0)
        fch_32(:,:,1)=0.0_R8
        fch_52(:,:,0)=0.0_R8
        fch_52(:,:,1)=fch(:,:,1)
      endif
      if(b2fstat_version.ge.'01.001.009') then
        call cfrure (nget, n2*ns, kinrgy, 'kinrgy')
      else
        do is = 0, ns-1
          kinrgy(:,:,is)=0.5_R8*am(is)*mp*
     .     (ua(:,:,is)**2+uadia(:,:,0,is)**2+uadia(:,:,1,is)**2)
        enddo
      endif
      call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
*   ..read time
      if(b2fstat_version.ge.'01.001.018') then                          !xpb
        call cfrure (nget, 1, time, 'time')
      else
        time(1)=0.0_R8
      endif
      tim=time(1)
      call xertst (0.0_R8.le.tim, 'faulty input time')
*   ..read parallel current
      if(b2fstat_version.ge.'01.001.069') then                          !srv 09.09.11
	   call cfrure (nget, n2, fch_p, 'fch_p')
	   fch_p(:,:,1)=0.0_R8                                            !srv 19.03.12
      elseif(b2fstat_version.ge.'01.001.025' .and.
     &   b2fstat_version.le.'01.001.068') then                          !xpb srv 10.06.08 {
        call cfrure (nget, n2*2, fch_p, 'fch_p')
	elseif (b2fstat_version.ge.'01.001.022' .and.
     &   b2fstat_version.ge.'01.001.024') then                          !srv 10.06.08 }
        call cfrure (nget, n2, fch_p, 'fch_p')
	   fch_p(:,:,1)=0.0_R8                                            !srv 19.03.12
      else
        fch_p(:,:,0) = fch(:,:,0)     ! this ensures the old behaviour is kept
        fch_p(:,:,1) = 0.0_R8         ! this ensures the old behaviour is kept
      endif                                                             !srv 10.06.08 }

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2rups

      end subroutine b2rups
