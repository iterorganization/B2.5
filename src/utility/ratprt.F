*-----------------------------------------------------------------------
*.specification

      subroutine ratprt (nprt, nz, npart, nidis, idis, njdis, jdis, 
     &  rtt, nt, rtn, nn, 
     &  rtsa, rtra, rtqa, rtqa_ra, rtcx, rtrd, rtbr, rtzz, rtz2, 
     &  rtpt, rtpi)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nz, npart, 
     *  nidis, idis(0:nidis-1), njdis, jdis(0:njdis-1),
     &  nt, nn
      real (kind=R8) ::
     *  rtt(0:nt), rtn(0:nn), rtsa(0:nt,0:nn,0:npart-1),
     &  rtra(0:nt,0:nn,0:npart-1), rtqa(0:nt,0:nn,0:npart-1),
     &  rtcx(0:nt,0:nn,0:npart-1), rtrd(0:nt,0:nn,0:npart-1),
     &  rtbr(0:nt,0:nn,0:npart-1), rtzz(0:nt,0:nn,0:npart-1),
     &  rtz2(0:nt,0:nn,0:npart-1), rtpt(0:nt,0:nn,0:npart-1),
     &  rtpi(0:nt,0:nn,0:npart-1), rtqa_ra(0:nt,0:nn,0:npart-1)
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     RATPRT prints an excerpt from a table of atomic rate
*     coefficients.
*
*
*  3. description
*
*     This subroutine writes to a formatted file an excerpt of a table
*     of atomic rate coefficients for ionisation, electron energy loss,
*     recombination, charge exchange, and radiation. The table provides
*     the rate coefficients as a function of density and temperature, but 
*     in this routine only one-dimensional data are displayed: the rate
*     coefficients as function of of temperature for several values
*     of the density, and as a function of density for several values
*     of the temperature.
*
*
*  5. parameters
*
*     nprt - integer, input.
*     nprt specifies a unit number for formatted output.
*
*     nz - integer, input.
*     nz specifies the nuclear charge of the species for which the rate
*     coefficients are supplied.
*
*     npart - integer, input.
*     npart specifies the number of partition of the homonuclear sequence.
*
*     nidis - integer, input.
*     nidis specifies the dimension of the array idis.
*
*     idis - (0:nidis-1) integer array.
*     idis is the selection of values of the temperature index at which
*     output will be produced to display the density dependence of the
*     rate coefficients.
*
*     njdis - integer, input.
*     njdis specifies the dimension of the array jdis.
*
*     jdis - (0:njdis-1) integer array.
*     jdis is the selection of values of the density index at which
*     output will be produced to display the temperature dependence of
*     the rate coefficients.
*
*     rtt - (0:nt) real*8 array, input.
*     rtt(0:nt) specifies the sequence of values of the temperature
*     variable, expressed in eV, for which the rate coefficients are
*     evaluated. rtt(0:nt) must be positive.
*
*     nt - integer, input.
*     nt specifies the upper dimension bound for the array rtt.
*     It will hold that 0.le.nt.
*
*     rtn - (0:nn) real*8 array, input.
*     rtn(0:nn) specifies the sequence of values of the density
*     variable, expressed in m**-3, for which the rate coefficients are
*     evaluated. rtn(0:nn) must be positive.
*
*     nn - integer, input.
*     nn specifies the upper dimension bound for the array rtn.
*     It will hold that 0.le.nn.
*
*     rtsa - (0:nt,0:nn,0:npart-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:npart-1), rtsa(i,j,k) specifies the
*     ionisation rate coefficient at electron temperature rtt(i) and
*     electron density rtn(j) for the process k->k+1. rtsa has
*     dimension m**3/s.
*
*     rtra - (0:nt,0:nn,0:npart-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:npart-1), rtra(i,j,k) specifies the
*     recombination rate coefficient at electron temperature rtt(i) and
*     electron density rtn(j) for the process k->k-1. rtra has
*     dimension m**3/s.
*
*     rtqa - (0:nt,0:nn,0:npart-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:npart-1), rtqa(i,j,k) specifies the
*     electron heat loss rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for processes starting from charge
*     state k. rtqa has dimension eV*m**3/s.
*
*     rtcx - (0:nt,0:nn,0:npart-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:npart-1), rtcx(i,j,k) specifies the charge
*     exchange rate coefficient at neutral hydrogen temperature rtt(i)
*     and electron density rtn(j) for the process k->k-1. rtcx has
*     dimension m**3/s.
*
*     rtrd - (0:nt,0:nn,0:npart-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:npart-1), rtrd(i,j,k) specifies the
*     line radiation rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for processes starting from charge
*     state k. rtrd has dimension eV*m**3/s.
*
*     rtbr - (0:nt,0:nn,0:npart-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:npart-1), rtbr(i,j,k) specifies the
*     bremsstrahlung radiation rate coefficient at electron temperature
*     rtt(i) and electron density rtn(j) emitted from charge
*     state k. rtbr has dimension eV*m**3/s.
*
*     rtzz - (0:nt,0:nn,0:npart-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:npart-1), rtzz(i,j,k) specifies the
*     effective charge at electron temperature rtt(i)
*     and electron density rtn(j) for stage k. 
*     rtzz has dimension e.
*
*     rtz2 - (0:nt,0:nn,0:npart-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:npart-1), rtz2(i,j,k) specifies the
*     effective square charge at electron temperature rtt(i)
*     and electron density rtn(j) for stage k. 
*     rtz2 has dimension e**2.
*
*     rtpt - (0:nt,0:nn,0:npart-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:npart-1), rtpt(i,j,k) specifies the
*     cumulative ionisation potential at electron temperature rtt(i)
*     and electron density rtn(j) to reach stage k. 
*     rtpt has dimension eV.
*
*     rtpi - (0:nt,0:nn,0:npart-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:npart-1), rtpi(i,j,k) specifies the
*     effective ionisation potential at electron temperature rtt(i)
*     and electron density rtn(j) to reach stage k from stage k-1. 
*     rtpi has dimension eV.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer i, j, k, ind
*   ..procedures
      external subini, subend, xertst

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
      call subini ('ratprt')
*   ..test nz, nt, nn
      call xertst (1.le.nz, 'ratprt--faulty argument nz')
      call xertst (1.le.npart.and.npart.le.nz+1,
     . 'ratprt--faulty argument npart')
      call xertst (0.le.nt.and.0.le.nn,
     . 'ratprt--faulty argument nt or nn')

* ..write summary of table
*   ..write basic data
      write (nprt,'(//2x,a,i2,a,i2,a)')
     &  '== nuclear charge number: ',nz,' with ',npart,' partitions =='
      write (nprt,101) 'rtt(0:nt)', (rtt(i),i=0,nt)
      write (nprt,101) 'rtn(0:nn)', (rtn(i),i=0,nn)
*   ..write rtsa as function of rtt for several rtn
      do k = 0, npart-2
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtsa(0:nt,',j,',',k,')',
     &    (rtsa(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtra as function of rtt for several rtn
      do k = 1, npart-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtra(0:nt,',j,',',k,')',
     &    (rtra(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtqa as function of rtt for several rtn
      do k = 0, npart-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtqa(0:nt,',j,',',k,')',
     &    (rtqa(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtqa_ra as function of rtt for several rtn
      do k = 0, npart-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtqa_ra(0:nt,',j,',',k,')',
     &    (rtqa_ra(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtcx as function of rtt for several rtn
      do k = 1, npart-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtcx(0:nt,',j,',',k,')',
     &    (rtcx(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtrd as function of rtt for several rtn
      do k = 0, npart-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtrd(0:nt,',j,',',k,')',
     &    (rtrd(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtbr as function of rtt for several rtn
      do k = 0, npart-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtbr(0:nt,',j,',',k,')',
     &    (rtbr(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtzz as function of rtt for several rtn
      do k = 0, npart-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtzz(0:nt,',j,',',k,')',
     &    (rtzz(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtz2 as function of rtt for several rtn
      do k = 0, npart-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtz2(0:nt,',j,',',k,')',
     &    (rtz2(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtpt as function of rtt for several rtn
      do k = 0, npart-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtpt(0:nt,',j,',',k,')',
     &    (rtpt(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtpi as function of rtt for several rtn
      do k = 0, npart-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtpi(0:nt,',j,',',k,')',
     &    (rtpi(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtsa as function of rtn for several rtt
      do k = 0, npart-2
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtsa(',i,',0:nn,',k,')',
     &    (rtsa(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtra as function of rtn for several rtt
      do k = 1, npart-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtra(',i,',0:nn,',k,')',
     &    (rtra(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtqa as function of rtn for several rtt
      do k = 0, npart-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtqa(',i,',0:nn,',k,')',
     &    (rtqa(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtcx as function of rtn for several rtt
      do k = 1, npart-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtcx(',i,',0:nn,',k,')',
     &    (rtcx(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtrd as function of rtn for several rtt
      do k = 0, npart-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtrd(',i,',0:nn,',k,')',
     &    (rtrd(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtbr as function of rtn for several rtt
      do k = 0, npart-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtbr(',i,',0:nn,',k,')',
     &    (rtbr(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtzz as function of rtn for several rtt
      do k = 0, npart-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtzz(',i,',0:nn,',k,')',
     &    (rtzz(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtz2 as function of rtn for several rtt
      do k = 0, npart-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtz2(',i,',0:nn,',k,')',
     &    (rtz2(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtpt as function of rtn for several rtt
      do k = 0, npart-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtpt(',i,',0:nn,',k,')',
     &    (rtpt(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtpi as function of rtn for several rtt
      do k = 0, npart-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtpi(',i,',0:nn,',k,')',
     &    (rtpi(i,j,k),j=0,nn)
       enddo
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.formats

  101 format (1p/2x,a/(12e11.3))
  102 format (1p/2x,a,i2,a,i2,a/(12e11.3))

*-----------------------------------------------------------------------
*.end ratprt

      end subroutine ratprt
