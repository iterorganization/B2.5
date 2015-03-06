*-----------------------------------------------------------------------
*.specification

      subroutine b2wups (nput, nx, ny, ns, 
     &  ne, na, ua, uadia, te, ti, po,
     &  fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy)
      use b2mod_types
      use b2mod_time
      use b2mod_version
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, nx, ny, ns
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  ne(-1:nx,-1:ny), uadia(-1:nx,-1:ny,0:1,0:ns-1),			 !xpb
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     1	fch_32(-1:nx,-1:ny,0:1), fch_52(-1:nx,-1:ny,0:1),
     2  fch_p(-1:nx,-1:ny,0:1), kinrgy(-1:nx,-1:ny,0:ns-1)         !srv 10.06.08
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WUPS writes the state of the plasma to an un*formatted file.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) the state of the plasma:
*     fields na, ua, te, ti, po, and fluxes fna, fhe, fhi, fch.
*     Unit nput must be connected to an un*formatted file. See the
*     code to discern the format of the output file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     nx, ny, ns - integers, input.
*     See the calling routine for description.
*
*     ne, na, ua, uadia, te, ti, po, fna, fhe, fhi, fch - real*8 arrays, input.
*     See the calling routine for description.
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2
      real (kind=R8) :: time(1)
*   ..procedures
      external subini, subend, xertst, cfwure
      external b2xvps, b2xvsg

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wups')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')

* ..save state and fluxes
*   ..save plasma state
      call cfwure (nput, n2*ns, na, 'na')
      if (newversion.ge.'01.001.024') 
     > call cfwure (nput, n2, ne, 'ne')
      call cfwure (nput, n2*ns, ua, 'ua')
      if (newversion.ge.'01.001.000')  
     > call cfwure (nput, n2*2*ns, uadia, 'uadia')
      call cfwure (nput, n2, te, 'te')
      call cfwure (nput, n2, ti, 'ti')
      call cfwure (nput, n2, po, 'po')
*   ..save fluxes
      call cfwure (nput, n2*2*ns, fna, 'fna')
      call cfwure (nput, n2*2, fhe, 'fhe')
      call cfwure (nput, n2*2, fhi, 'fhi')
      call cfwure (nput, n2*2, fch, 'fch')
      if (newversion.ge.'01.001.005')  
     > call cfwure (nput, n2*2, fch_32, 'fch_32')
      if (newversion.ge.'01.001.005')  
     > call cfwure (nput, n2*2, fch_52, 'fch_52')
      if (newversion.ge.'01.001.009')  
     > call cfwure (nput, n2*ns, kinrgy, 'kinrgy')
*   ..save time
      time(1)=tim
      if (newversion.ge.'01.001.018')  
     > call cfwure (nput, 1, time, 'time')
*   ..save parallel current
      if(newversion.ge.'01.001.069') then                       !srv 09.09.11
	 call cfwure (nput, n2, fch_p(:,:,0), 'fch_p')
      elseif(newversion.ge.'01.001.025' .and.
     &   newversion.le.'01.001.068') then                       !srv 10.06.08 {
        call cfwure (nput, n2*2, fch_p, 'fch_p')
      elseif (newversion.ge.'01.001.022' .and.
     &   newversion.ge.'01.001.024') then
        call cfwure (nput, n2, fch_p(:,:,0), 'fch_p')
      endif                                                     !srv 10.06.08 }
* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wups

      end
