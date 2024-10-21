*-----------------------------------------------------------------------
*.specification

      subroutine b2wfrc (nprt, ns0, nidis, idis, njdis, jdis, rtt, nt,
     &  rtn, nn, zamin, zamax, zn, fix_recomb, zmax_recomb,
     &  rtlsa, rtlra, rtlqa, rtlqr, rtlcx, rtlza, rtlz2, rtlpt,
     &  rtlpi)
      use b2mod_types
      use b2mod_math
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, ns0, nidis, idis(0:nidis-1), njdis, 
     &  jdis(0:njdis-1), nt, nn, fix_recomb, zmax_recomb
      real (kind=R8) ::
     *  rtt(0:nt), rtn(0:nn), 
     &  zamin(0:ns0-1), zamax(0:ns0-1), zn(0:ns0-1),
     &  rtlsa(0:nt,0:nn,0:ns0-1), rtlra(0:nt,0:nn,0:ns0-1),
     &  rtlqa(0:nt,0:nn,0:ns0-1), rtlcx(0:nt,0:nn,0:ns0-1),
     &  rtlza(0:nt,0:nn,0:ns0-1), rtlz2(0:nt,0:nn,0:ns0-1),
     &  rtlpt(0:nt,0:nn,0:ns0-1), rtlpi(0:nt,0:nn,0:ns0-1),
     &  rtlqr(0:nt,0:nn,0:ns0-1)
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WFRC writes to a formatted file a table of atomic rate
*     coefficients that has been prepared for B2.5.
*
*
*  3. description
*
*     This subroutine writes to unit nprt a table of atomic rate
*     coefficients for ionisation, electron energy loss, recombination,
*     and charge exchange. Unit nprt must be connected to a formatted
*     file.
*
*     The table has a two-dimensional structure, with density and
*     temperature as independent variables, but in this routine only
*     one-dimensional data are displayed: the rate coefficients as
*     function of of temperature for several values of the density,
*     and as a function of density for several values of the
*     temperature.
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
*     zamin, zamax - (0:ns0-1) real*8 arrays, input.
*     For k in (0:ns0-1), zamin(k) and zamax(k) specify the atomic charge 
*     range of species k. zamin and zamax are dimensionless.
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
*     charge state k, excluding the contributions 
*     coming from potential energy of recombination and other recombination
*     process contributions. In short, it includes the processes not
*     accounted for in Eirene. rtqa has dimension eV*m**3/s.
*
*     rtlqr - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), expu(rtlqr(i,j,k)) specifies
*     the electron heat loss rate coefficient at electron temperature
*     rtt(i) and electron density rtn(j) for processes starting from
*     charge state k. rtqr only includes the contribution of recombination
*     which replaces AMJUEL when Eirene is not used.
*     rtqr has dimension eV*m**3/s.
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
*     and electron density rtn(j) for the superstage k. rtlza has
*     dimension e.
*
*     rtlz2 - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), rtlz2(i,j,k) specifies
*     the effective square charge at electron temperature rtt(i)
*     and electron density rtn(j) for the superstage k. rtlz2 has
*     dimension e**2.
*
*     rtlpt - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), rtlpt(i,j,k) specifies
*     the cumulative ionisation potential at electron temperature rtt(i)
*     and electron density rtn(j) to reach superstage k. rtlpt has
*     dimension eV.
*
*     rtlpi - (0:nt,0:nn,0:ns0-1) real*8 array, input.
*     For (i,j,k) in (0:nt,0:nn,0:ns0-1), rtlpi(i,j,k) specifies
*     the effective ionisation potential at electron temperature rtt(i)
*     and electron density rtn(j) to reach superstage k from
*     superstage k-1. rtlpi has dimension eV.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer i, j, k, ind
*   ..procedures
      external xertst, b2xvsg

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wfrc')
*   ..test ns0, nt, nn
      call xertst (1.le.ns0, 'faulty argument ns0')
      call xertst (0.le.nt.and.0.le.nn, 'faulty argument nt or nn')
*   ..test sign of rtt, rtn
      call b2xvsg (nt+1, rtt, 1, 'rtt', '.gt.')
      call b2xvsg (nn+1, rtn, 1, 'rtn', '.gt.')
*   ..set cutll
      if (.not.b2mod_math_initialised) call b2mod_math_init

* ..write summary of table
*   ..write rtsa as function of rtt for several rtn
      do k = 0, ns0-2
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtsa(0:nt,',j,',',k,')',
     &    zamin(k), zamax(k), zn(k), (expu1(rtlsa(i,j,k)),i=0,nt)
       enddo
      enddo
*   ..write rtsa as function of rtn for several rtt
      do k = 0, ns0-2
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtsa(',i,',0:nn,',k,')',
     &    zamin(k), zamax(k), zn(k), (expu1(rtlsa(i,j,k)),j=0,nn)
       enddo
      enddo
*   ..write rtra as function of rtt for several rtn
      do k = 1, ns0-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtra(0:nt,',j,',',k,')',
     &    zamin(k), zamax(k), zn(k), (expu1(rtlra(i,j,k)),i=0,nt)
       enddo
      enddo
*   ..write rtra as function of rtn for several rtt
      do k = 1, ns0-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtra(',i,',0:nn,',k,')',
     &    zamin(k), zamax(k), zn(k), (expu1(rtlra(i,j,k)),j=0,nn)
       enddo
      enddo
*   ..write rtqa as function of rtt for several rtn
      do k = 0, ns0-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtqa(0:nt,',j,',',k,')',
     &    zamin(k), zamax(k), zn(k), (expu1(rtlqa(i,j,k)),i=0,nt)
       enddo
      enddo
*   ..write rtqa as function of rtn for several rtt
      do k = 0, ns0-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtqa(',i,',0:nn,',k,')',
     &    zamin(k), zamax(k), zn(k), (expu1(rtlqa(i,j,k)),j=0,nn)
       enddo
      enddo
*   ..write rtqa (with fix_recomb) as function of rtn for several rtt
      if (fix_recomb.ne.0) then
       do k = 0, ns0-1
        if (nint(zn(k)).gt.zmax_recomb) cycle
        do ind = 0, nidis-1
         i = idis(ind)
         write (nprt,102) 'rtqa[with fix_recomb](',i,',0:nn,',k,')',
     &    zamin(k), zamax(k), zn(k),
     &   (expu1(rtlqa(i,j,k))-expu1(rtlra(i,j,k))*rtlpi(i,j,k),j=0,nn)
        enddo
       enddo
      endif
*   ..write rtqr as function of rtt for several rtn
      do k = 0, ns0-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtqr(0:nt,',j,',',k,')',
     &    zamin(k), zamax(k), zn(k), (expu1(rtlqr(i,j,k)),i=0,nt)
       enddo
      enddo
*   ..write rtqr as function of rtn for several rtt
      do k = 0, ns0-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtqr(',i,',0:nn,',k,')',
     &    zamin(k), zamax(k), zn(k), (expu1(rtlqr(i,j,k)),j=0,nn)
       enddo
      enddo
*   ..write rtcx as function of rtt for several rtn
      do k = 1, ns0-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtcx(0:nt,',j,',',k,')',
     &    zamin(k), zamax(k), zn(k), (expu1(rtlcx(i,j,k)),i=0,nt)
       enddo
      enddo
*   ..write rtcx as function of rtn for several rtt
      do k = 1, ns0-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtcx(',i,',0:nn,',k,')',
     &    zamin(k), zamax(k), zn(k), (expu1(rtlcx(i,j,k)),j=0,nn)
       enddo
      enddo
*   ..write rtza as function of rtt for several rtn
      do k = 0, ns0-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtza(0:nt,',j,',',k,')',
     &    zamin(k), zamax(k), zn(k), (rtlza(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtza as function of rtn for several rtt
      do k = 0, ns0-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtza(',i,',0:nn,',k,')',
     &    zamin(k), zamax(k), zn(k), (rtlza(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtz2 as function of rtt for several rtn
      do k = 0, ns0-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtz2(0:nt,',j,',',k,')',
     &    zamin(k), zamax(k), zn(k), (rtlz2(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtz2 as function of rtn for several rtt
      do k = 0, ns0-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtz2(',i,',0:nn,',k,')',
     &    zamin(k), zamax(k), zn(k), (rtlz2(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtpt as function of rtt for several rtn
      do k = 0, ns0-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtpt(0:nt,',j,',',k,')',
     &    zamin(k), zamax(k), zn(k), (rtlpt(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtpt as function of rtn for several rtt
      do k = 0, ns0-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtpt(',i,',0:nn,',k,')',
     &    zamin(k), zamax(k), zn(k), (rtlpt(i,j,k),j=0,nn)
       enddo
      enddo
*   ..write rtpi as function of rtt for several rtn
      do k = 0, ns0-1
       do ind = 0, njdis-1
        j = jdis(ind)
        write (nprt,102) 'rtpi(0:nt,',j,',',k,')',
     &    zamin(k), zamax(k), zn(k), (rtlpi(i,j,k),i=0,nt)
       enddo
      enddo
*   ..write rtpi as function of rtn for several rtt
      do k = 0, ns0-1
       do ind = 0, nidis-1
        i = idis(ind)
        write (nprt,102) 'rtpi(',i,',0:nn,',k,')',
     &    zamin(k), zamax(k), zn(k), (rtlpi(i,j,k),j=0,nn)
       enddo
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.formats

  102 format (0p/2x,a,i3,a,i3,a,5x,
     &  '(zamin =',f6.2,3x,'zamax =',f6.2,3x,'zn =',f6.2,')'/
     &  (1p,12e11.3))

*-----------------------------------------------------------------------
*.end b2wfrc

      end subroutine b2wfrc
