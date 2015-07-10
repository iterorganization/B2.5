*-----------------------------------------------------------------------
*.specification

      subroutine b2wgrc (nprt, ns0, nidis, idis, njdis, jdis, rtt, nt,
     &  rtn, nn, zn, 
     &  rtlsa, rtlra, rtlqa, rtlcx, rtlza, rtlz2, rtlpt, rtlpi)
      use b2mod_types
      use b2mod_constants
      implicit none
#     include "emissivity.constants"
*   ..input arguments (unchanged on exit)
      integer nprt, ns0, nidis, idis(0:nidis-1), njdis, jdis(0:njdis-1),
     &  nt, nn
      real (kind=R8) ::
     *   rtt(0:nt), rtn(0:nn), zn(0:ns0-1),
     &  rtlsa(0:nt,0:nn,0:ns0-1), rtlra(0:nt,0:nn,0:ns0-1),
     &  rtlqa(0:nt,0:nn,0:ns0-1), rtlcx(0:nt,0:nn,0:ns0-1),
     &  rtlza(0:nt,0:nn,0:ns0-1), rtlz2(0:nt,0:nn,0:ns0-1),
     &  rtlpt(0:nt,0:nn,0:ns0-1), rtlpi(0:nt,0:nn,0:ns0-1)
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGRC displays a table of atomic rate coefficients that has been
*     prepared for b2.
*
*
*  3. description
*
*     This subroutine displays a table of atomic rate coefficients
*     for ionisation, electron energy loss, recombination, and
*     charge-exchange. The table has a two-dimensional structure,
*     with density and temperature as independent variables, but
*     in this routine only one-dimensional data are displayed:
*     the rate coefficients as function of of temperature for several
*     values of the density, and as a function of density for several
*     values of the temperature.
*
*
*  5. parameters
*
*     nprt - integer, input.
*     nprt specifies a unit number for formatted output.
*
*     ns0 - integer, input.
*     ns0 specifies the number of species. It will hold that 1.le.ns0.
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
*     variable, expressed in ev, for which the rate coefficients are
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
*     zn - (0:ns0-1) real*8 array, input.
*     For k in (0:ns0-1), zn(k) specifies the nuclear charge of species
*     k. zn is dimensionless.
*
*     rtlsa - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), expu(rtlsa(i,j,k)) specifies
*     the ionisation rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for the process k->k+1. rtsa has
*     dimension m**3/s.
*
*     rtrla - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), expu(rtlra(i,j,k)) specifies
*     the recombination rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for the process k->k-1. rtra has
*     dimension m**3/s.
*
*     rtlqa - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), expu(rtlqa(i,j,k)) specifies
*     the electron heat loss rate coefficient at electron temperature
*     rtt(i) and electron density rtn(j) for processes starting from
*     charge state k. rtqa has dimension ev*m**3/s.
*
*     rtlcx - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), expu(rtlcx(i,j,k)) specifies
*     the charge exchange rate coefficient at neutral hydrogen
*     temperature rtt(i) and electron density rtn(j) for the process
*     k->k-1. rtcx has dimension m**3/s.
*
*     rtlza - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), rtlza(i,j,k) specifies
*     the effective charge state at electron temperature rtt(i)
*     and electron density rtn(j) for the superstage k. rtza has
*     dimension e.
*
*     rtlz2 - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), rtlz2(i,j,k) specifies
*     the effective square charge at electron temperature rtt(i)
*     and electron density rtn(j) for the superstage k. rtz2 has
*     dimension e**2.
*
*     rtlpt - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), rtlpt(i,j,k) specifies
*     the cumulative ionisation potential at electron temperature rtt(i)
*     and electron density rtn(j) to reach superstage k. rtpt has
*     dimension eV.
*
*     rtlpi - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), rtlpi(i,j,k) specifies
*     the effective ionisation potential at electron temperature rtt(i)
*     and electron density rtn(j) to reach superstage k from
*     superstage k-1. rtpi has dimension eV.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nmax, i, j, k, ind, izn0
      character header*60, chk*2, chdc*(*), fmtp*(*)
      parameter (chdc='''$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$''',
     &  fmtp='(2x,a6,4x,a)')
      real (kind=R8) ::
     *   cutll, rydberg, Ei, Ef, t0
      real (kind=R8), allocatable :: wt(:), wn(:)
*   ..procedures
      intrinsic nint, max, exp, log
      real (kind=R8) ::
     *   machsfr, expu1
      real (kind=R4) small_real_constant
      external subini, subend, xertst, pagini, pagend, machsfr,
     &  dashdc, curved_d
      expu1(t0) = exp(max(t0,cutll))+small_real_constant
*     (expu1 is the exponential function with underflow set to a tiny
*     positive value.)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgrc')
*   ..test nprt
      call xertst (1.le.nprt, 'faulty argument nprt')
*   ..test ns0, nt, nn
      call xertst (1.le.ns0, 'faulty argument ns0')
      call xertst (0.le.nt.and.0.le.nn, 'faulty argument nt or nn')
      allocate (wn(0:nn))
      allocate (wt(0:nt))
*   ..set cutll
      cutll = log(2*machsfr())
      small_real_constant=tiny(small_real_constant)
*!!   (increase cutll if a floating point error occurs.)
      rydberg=(me/(2.0_R8*qe))*(qe**2/(2.0_R8*const_h*eps0))**2
*   ..set izn0
      izn0 = nint(zn(0))
*   ..compute ionisation potential scale
      Ef = 1.10_R8 * zn(0)**2
      Ei = 1.12_R8 * zn(0)**2.4_R8
*   ..write blank line
      write (nprt,'()')

* ..plot summary of table
*!!   Calls to dashdc have been commented out; they give rise to a
*     floating point exception on the Sparc. An NCARG error ?
*   ..plot rtsa as function of rtt for several rtn
      do ind = 0, njdis-1
       j = jdis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'IONISATION RATE COEFS,        Z0 =',izn0,', NE =',rtn(j)
       write (nprt,fmtp) '+rtsa', header
       call pagini (3, rtt(0), rtt(nt), 1.0e-20_R8, 1.0e-12_R8, 'rtt',
     &   'rtsa', header)
       do k = 0, ns0-2
        do i = 0, nt
         wt(i) = expu1(rtlsa(i,j,k))
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	if(real(maxval(wt(0:nt)),R4).gt.0) then
	  call curved_d (rtt, wt, nt+1)
	endif
       enddo
       call pagend ()
      enddo
*   ..plot rtsa as function of rtn for several rtt
      do ind = 0, nidis-1
       i = idis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'IONISATION RATE COEFS,        Z0 =',izn0,', TE =',rtt(i)
       write (nprt,fmtp) '+rtsa', header
       call pagini (3, rtn(0), rtn(nn), 1.0e-20_R8, 1.0e-12_R8, 'rtn',
     &   'rtsa', header)
       do k = 0, ns0-2
        do j = 0, nn
         wn(j) = expu1(rtlsa(i,j,k))
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	if(real(maxval(wn(0:nn)),R4).gt.0) then
	  call curved_d (rtn, wn, nn+1)
	endif
       enddo
       call pagend ()
      enddo
*   ..plot rtra as function of rtt for several rtn
      do ind = 0, njdis-1
       j = jdis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'RECOMBINATION RATE COEFS,     Z0 =',izn0,', NE =',rtn(j)
       write (nprt,fmtp) '+rtra', header
       call pagini (3, rtt(0), rtt(nt), 1.0e-22_R8, 1.0e-14_R8, 'rtt',
     &   'rtra', header)
       do k = 1, ns0-1
        do i = 0, nt
         wt(i) = expu1(rtlra(i,j,k))
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	if(real(maxval(wt(0:nt)),R4).gt.0) then
	  call curved_d (rtt, wt, nt+1)
	endif
       enddo
       call pagend ()
      enddo
*   ..plot rtra as function of rtn for several rtt
      do ind = 0, nidis-1
       i = idis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'RECOMBINATION RATE COEFS,     Z0 =',izn0,', TE =',rtt(i)
       write (nprt,fmtp) '+rtra', header
       call pagini (3, rtn(0), rtn(nn), 1.0e-22_R8, 1.0e-14_R8, 'rtn',
     &   'rtra', header)
       do k = 1, ns0-1
        do j = 0, nn
         wn(j) = expu1(rtlra(i,j,k))
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	if(real(maxval(wn(0:nn)),R4).gt.0) then
	  call curved_d (rtn, wn, nn+1)
	endif
       enddo
       call pagend ()
      enddo
*   ..plot rtqa as function of rtt for several rtn
      do ind = 0, njdis-1
       j = jdis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'ENERGY LOSS RATE COEFS,       Z0 =',izn0,', NE =',rtn(j)
       write (nprt,fmtp) '+rtqa', header
       call pagini (3, rtt(0), rtt(nt), 1.0e-18_R8, 1.0e-10_R8, 'rtt',
     &   'rtqa', header)
       do k = 0, ns0-1
        do i = 0, nt
         wt(i) = expu1(rtlqa(i,j,k))
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	if(real(maxval(wt(0:nt)),R4).gt.0) then
	  call curved_d (rtt, wt, nt+1)
	endif
       enddo
       call pagend ()
      enddo
*   ..plot rtqa as function of rtn for several rtt
      do ind = 0, nidis-1
       i = idis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'ENERGY LOSS RATE COEFS,       Z0 =',izn0,', TE =',rtt(i)
       write (nprt,fmtp) '+rtqa', header
       call pagini (3, rtn(0), rtn(nn), 1.0e-18_R8, 1.0e-10_R8, 'rtn',
     &   'rtqa', header)
       do k = 0, ns0-1
        do j = 0, nn
         wn(j) = expu1(rtlqa(i,j,k))
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	if(real(maxval(wn(0:nn)),R4).gt.0) then
	  call curved_d (rtn, wn, nn+1)
	endif
       enddo
       call pagend ()
      enddo
*   ..plot rtcx as function of rtt for several rtn
      do ind = 0, njdis-1
       j = jdis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'CX-RECOMBINATION RATE COEFS,  Z0 =',izn0,', NN =',rtn(j)
       write (nprt,fmtp) '+rtcx', header
       call pagini (3, rtt(0), rtt(nt), 1.0e-18_R8, 1.0e-10_R8, 'rtt',
     &   'rtcx', header)
       do k = 1, ns0-1
        do i = 0, nt
         wt(i) = expu1(rtlcx(i,j,k))
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	if(real(maxval(wt(0:nt)),R4).gt.0) then
	  call curved_d (rtt, wt, nt+1)
	endif
       enddo
       call pagend ()
      enddo
*   ..plot rtcx as function of rtn for several rtt
      do ind = 0, nidis-1
       i = idis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'CX-RECOMBINATION RATE COEFS,  Z0 =',izn0,', TI =',rtt(i)
       write (nprt,fmtp) '+rtcx', header
       call pagini (3, rtn(0), rtn(nn), 1.0e-18_R8, 1.0e-10_R8, 'rtn',
     &   'rtcx', header)
       do k = 1, ns0-1
        do j = 0, nn
         wn(j) = expu1(rtlcx(i,j,k))
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	if(real(maxval(wn(0:nn)),R4).gt.0) then
	  call curved_d (rtn, wn, nn+1)
	endif
       enddo
       call pagend ()
      enddo
*   ..plot rtza as function of rtt for several rtn
      do ind = 0, njdis-1
       j = jdis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'EFFECTIVE CHARGE,             Z0 =',izn0,', NE =',rtn(j)
       write (nprt,fmtp) '+rtza', header
       call pagini (6, rtt(0), rtt(nt), -0.5_R8, zn(0)+0.5_R8, 'rtt',
     &   'rtza', header)
       do k = 0, ns0-1
        do i = 0, nt
         wt(i) = rtlza(i,j,k)
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	call curved_d (rtt, wt, nt+1)
       enddo
       call pagend ()
      enddo
*   ..plot rtza as function of rtn for several rtt
      do ind = 0, nidis-1
       i = idis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'EFFECTIVE CHARGE,             Z0 =',izn0,', TE =',rtt(i)
       write (nprt,fmtp) '+rtza', header
       call pagini (6, rtn(0), rtn(nn), -0.5_R8, zn(0)+0.5_R8, 'rtn',
     &   'rtza', header)
       do k = 0, ns0-1
        do j = 0, nn
         wn(j) = rtlza(i,j,k)
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	call curved_d (rtn, wn, nn+1)
       enddo
       call pagend ()
      enddo
*   ..plot rtz2 as function of rtt for several rtn
      do ind = 0, njdis-1
       j = jdis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'EFFECTIVE SQUARE CHARGE,      Z0 =',izn0,', NE =',rtn(j)
       write (nprt,fmtp) '+rtz2', header
       call pagini (6, rtt(0), rtt(nt), -0.5_R8, zn(0)**2+0.5_R8, 'rtt',
     &   'rtz2', header)
       do k = 0, ns0-1
        do i = 0, nt
         wt(i) = rtlz2(i,j,k)
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	call curved_d (rtt, wt, nt+1)
       enddo
       call pagend ()
      enddo
*   ..plot rtz2 as function of rtn for several rtt
      do ind = 0, nidis-1
       i = idis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'EFFECTIVE SQUARE CHARGE,      Z0 =',izn0,', TE =',rtt(i)
       write (nprt,fmtp) '+rtz2', header
       call pagini (6, rtn(0), rtn(nn), -0.5_R8, zn(0)**2+0.5_R8, 'rtn',
     &   'rtz2', header)
       do k = 0, ns0-1
        do j = 0, nn
         wn(j) = rtlz2(i,j,k)
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	call curved_d (rtn, wn, nn+1)
       enddo
       call pagend ()
      enddo
*   ..plot rtpt as function of rtt for several rtn
      do ind = 0, njdis-1
       j = jdis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'CUMULATIVE IONISATION POT.  , Z0 =',izn0,', NE =',rtn(j)
       write (nprt,fmtp) '+rtpt', header
       call pagini (6, rtt(0), rtt(nt), -0.5_R8, Ei*rydberg+0.5_R8, 
     &   'rtt', 'rtpt', header)
       do k = 0, ns0-1
        do i = 0, nt
         wt(i) = rtlpt(i,j,k)
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	call curved_d (rtt, wt, nt+1)
       enddo
       call pagend ()
      enddo
*   ..plot rtpt as function of rtn for several rtt
      do ind = 0, nidis-1
       i = idis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'CUMULATIVE IONISATION POT.  , Z0 =',izn0,', TE =',rtt(i)
       write (nprt,fmtp) '+rtpt', header
       call pagini (6, rtn(0), rtn(nn), -0.5_R8, Ei*rydberg+0.5_R8, 
     &   'rtn', 'rtpt', header)
       do k = 0, ns0-1
        do j = 0, nn
         wn(j) = rtlpt(i,j,k)
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	call curved_d (rtn, wn, nn+1)
       enddo
       call pagend ()
      enddo
*   ..plot rtpi as function of rtt for several rtn
      do ind = 0, njdis-1
       j = jdis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'EFFECTIVE IONISATION POT.,    Z0 =',izn0,', NE =',rtn(j)
       write (nprt,fmtp) '+rtpi', header
       call pagini (6, rtt(0), rtt(nt), -0.5_R8, Ef*rydberg+0.5_R8, 
     &   'rtt', 'rtpi', header)
       do k = 0, ns0-1
        do i = 0, nt
         wt(i) = rtlpi(i,j,k)
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	call curved_d (rtt, wt, nt+1)
       enddo
       call pagend ()
      enddo
*   ..plot rtpi as function of rtn for several rtt
      do ind = 0, nidis-1
       i = idis(ind)
       write (header,'(a,i2,a,1pe9.2)')
     &   'EFFECTIVE IONISATION POT.,    Z0 =',izn0,', TE =',rtt(i)
       write (nprt,fmtp) '+rtpi', header
       call pagini (6, rtn(0), rtn(nn), -0.5_R8, Ef*rydberg+0.5_R8, 
     &   'rtn', 'rtpi', header)
       do k = 0, ns0-1
        do j = 0, nn
         wn(j) = rtlpi(i,j,k)
        enddo
        write (chk,'(i2)') k
*$$$    call dashdc (chdc//chk, 10, 1)
	call curved_d (rtn, wn, nn+1)
       enddo
       call pagend ()
      enddo

* ..return
      deallocate(wn,wt)
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgrc

      end subroutine b2wgrc
