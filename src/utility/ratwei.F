*-----------------------------------------------------------------------
*.specification

      subroutine ratwei (nrid, nz, rtt, nt, rtn, nn, rtsa, rtra, rtqa,
     &  rtcx, rtrd, rtbr)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nrid, nz, nt, nn
      real (kind=R8) ::
     *  rtt(0:nt), rtn(0:nn)
*   ..input/output arguments
      real (kind=R8) ::
     *  rtsa(0:nt,0:nn,0:nz), rtra(0:nt,0:nn,0:nz),
     &  rtqa(0:nt,0:nn,0:nz), rtcx(0:nt,0:nn,0:nz),
     &  rtrd(0:nt,0:nn,0:nz), rtbr(0:nt,0:nn,0:nz)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     RATWEI uses WEISHEIT data to improve tables of rate coefficients
*     for ionisation, recombination, electron energy loss, charge
*     exchange, and radiation; as function of temperature and density.
*     The WEISHEIT data set is incomplete.  The call to ratwei must
*     be preceded by a call to one of the comprehensive routines.
*
*
*  3. description
*
*     This routine employs WEISHEIT data /1,2/ to improve a table
*     of atomic rate coefficients.  These data were made available
*     by D. K. Stotler and D. E. Post of Princeton PPL.
*
*     rtt represents temperature or temperature per amu; units are ev.
*     rtn represents density; units are m**(-3).
*     rtsa(i,j,k) is ionisation rate coefficient; m**3/s.
*     rtra(i,j,k) is recombination rate coefficient; m**3/s.
*     rtqa(i,j,k) is energy loss rate coefficient; eV*m**3/s.
*     rtcx(i,j,k) is charge exchange rate coefficient; m**3/s.
*     rtrd(i,j,k) is line radiation rate coefficient; eV*m**3/s.
*     rtbr(i,j,k) is bremsstrahlung radiation rate coefficient; eV*m**3/s.
*     all at coordinates rtt(i) and rtn(j), and for the transition
*     starting from charge state k.
*     For rtsa, rtra, rtqa, rtrd, and rtbr, the quantities rtt and rtn are
*     interpreted as the electron temperature and density. For rtcx,
*     the quantity rtt is interpreted as the neutral hydrogen
*     temperature per atomic mass, and rtn is interpreted as the
*     electron density.
*
*     The Weisheit table provides data for hydrogen only, and only for
*     the processes of ionisation, recombination, and neutral hydrogen
*     energy loss and radiation.  The rate coefficients for charge exchange 
*     and for charged hydrogen energy loss must be obtained by a preceding 
*     call to one of the comprehensive packages, such as STRAHL or ADPAK.
*     The call to the present routine will then cause the rate
*     coefficients for the selected processes to be overwritten.
*
*  4. references
*
*     /1/ Jon C. Weisheit, "Recombination in Dense Plasmas", J. Phys. B
*     (Atom. Molec. Phys.) 8 (1975) 2556.
*     /2/ R. K. Janev, D. E. Post, W. D. Langer, K. Evans, D. B.
*     Heifetz and J. C. Weisheit, "Survey of Atomic Processes in Edge
*     Plasmas", J. Nuclear Mater. 121 (1984) 10-16.
*
*
*  5. parameters
*
*     nrid - integer, input.
*     nrid specifies a unit number for formatted input.
*
*     nz - integer, input.
*     nz specifies the nuclear charge of the species for which the rate
*     coefficients are to be evaluated.  In the present code it is
*     required that nz.eq.1.
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
*     rtsa - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtsa(i,j,k) specifies the
*     ionisation rate coefficient at electron temperature rtt(i) and
*     electron density rtn(j) for the process k->k+1. rtsa has
*     dimension m**3/s.
*     This subroutine will override rtsa(,,0).
*
*     rtra - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtra(i,j,k) specifies the
*     recombination rate coefficient at electron temperature rtt(i) and
*     electron density rtn(j) for the process k->k-1. rtra has
*     dimension m**3/s.
*     This subroutine will override rtra(,,nz).
*
*     rtqa - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtqa(i,j,k) specifies the
*     electron heat loss rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for processes starting from charge
*     state k. rtqa has dimension eV*m**3/s.
*     This subroutine will override rtqa(,,0).
*
*     rtcx - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtcx(i,j,k) specifies the charge
*     exchange rate coefficient at neutral hydrogen temperature rtt(i)
*     and electron density rtn(j) for the process k->k-1. rtcx has
*     dimension m**3/s.
*     This subroutine will leave rtcx unchanged.
*
*     rtrd - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtrd(i,j,k) specifies the
*     line radiation rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for processes starting from charge
*     state k. rtrd has dimension eV*m**3/s.
*     This subroutine will override rtrd(,,0).
*
*     rtbr - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtbr(i,j,k) specifies the
*     bremsstrahlung radiation rate coefficient at electron temperature 
*     rtt(i) and electron density rtn(j) emitted from charge
*     state k. rtbr has dimension eV*m**3/s.
*     This subroutine will leave rtbr(,,0) unchanged.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer mpd, mpe, i, j, iw, jw, n3
      parameter (mpd=15, mpe=60)
      real (kind=R8) ::
     *  wsveh(mpe,mpd), wsveh0(mpe,mpd), welms(mpe,mpd)
      real (kind=R8) ::
     *  cm3, fxt, fxn, t0, t1
      parameter (cm3=1.0e-6_R8)
*   ..procedures
      intrinsic exp, log
      external subini, subend, xertst, b2xvsg
      real (kind=R8) ::
     *   dens, temp, wsven
      dens(jw) = 10.0_R8**(16.0_R8+0.5_R8*(jw-1))
      temp(iw) = 10.0_R8**(0.1_R8*(iw-13))
      wsven(iw,jw) = wsveh(iw,jw)*welms(iw,jw)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
      call subini ('ratwei')
*   ..set n3
      n3 = (nz+1)*(nt+1)*(nn+1)
*   ..test nz, nt, nn
      call xertst (nz.eq.1, 'faulty argument nz')
      call xertst (0.le.nt.and.0.le.nn, 'faulty argument nt, nn')
*   ..test rtt, rtn
      call b2xvsg (nt+1, rtt, 1, 'rtt', '.gt.')
      call b2xvsg (nn+1, rtn, 1, 'rtn', '.gt.')
*   ..test rtsa, rtra, rtqa, rtcx, rtrd, rtbr
      call b2xvsg (n3, rtsa, 1, 'rtsa', '.ge.')
      call b2xvsg (n3, rtra, 1, 'rtra', '.ge.')
      call b2xvsg (n3, rtqa, 1, 'rtqa', '.ge.')
      call b2xvsg (n3, rtcx, 1, 'rtcx', '.ge.')
      call b2xvsg (n3, rtrd, 1, 'rtrd', '.ge.')
      call b2xvsg (n3, rtbr, 1, 'rtbr', '.ge.')
*   ..read data
      do jw = 1, mpd
       read (nrid,'(10(6(1x,e12.5)/))') (wsveh(iw,jw),iw=1,mpe)
      enddo
      do jw = 1, mpd
       read (nrid,'(10(6(1x,e12.5)/))') (wsveh0(iw,jw),iw=1,mpe)
      enddo
      do jw = 1, mpd
       read (nrid,'(10(6(1x,e12.5)/))') (welms(iw,jw),iw=1,mpe)
      enddo

* ..correct selected coefficients
*   ..initialise iw (temperature point in Weisheit table)
      iw = 1
*   ..loop over temperature
      do i = 0, nt
*    ..update iw
    1  if (iw.lt.mpe-1) then
        if (temp(iw+1).lt.rtt(i)) then
         iw = iw+1
         goto 1
        endif
       endif
*    ..compute coefficient for interpolation or extrapolation
       fxt = (log(rtt(i))-log(temp(iw)))/
     &   (log(temp(iw+1))-log(temp(iw)))
*    ..initialise jw (density point in Weisheit table)
       jw = 1
*    ..loop over density
       do j = 0, nn
*     ..update jw
    2   if (jw.lt.mpd-1) then
         if (dens(jw+1).lt.rtn(j)) then
          jw = jw+1
          goto 2
         endif
        endif
*     ..compute coefficient for interpolation or extrapolation
        fxn = (log(rtn(j))-log(dens(jw)))/
     &    (log(dens(jw+1))-log(dens(jw)))
*     ..compute new rate coefficients
        t0 = (1.0_R8-fxt)*log(wsveh(iw,jw))+fxt*log(wsveh(iw+1,jw))
        t1 = (1.0_R8-fxt)*log(wsveh(iw,jw+1))+fxt*log(wsveh(iw+1,jw+1))
        rtsa(i,j,0) = cm3*exp((1.0_R8-fxn)*t0+fxn*t1)
        t0 = (1.0_R8-fxt)*log(wsveh0(iw,jw))+fxt*log(wsveh0(iw+1,jw))
        t1 = (1.0_R8-fxt)*log(wsveh0(iw,jw+1))+
     .                fxt*log(wsveh0(iw+1,jw+1))
        rtra(i,j,nz) = cm3*exp((1.0_R8-fxn)*t0+fxn*t1)
        t0 = (1.0_R8-fxt)*log(wsven(iw,jw))+fxt*log(wsven(iw+1,jw))
        t1 = (1.0_R8-fxt)*log(wsven(iw,jw+1))+fxt*log(wsven(iw+1,jw+1))
        rtrd(i,j,0) = rtrd(i,j,0)*
     .    (cm3*exp((1.0_R8-fxn)*t0+fxn*t1))/rtqa(i,j,0)
        rtqa(i,j,0) = cm3*exp((1.0_R8-fxn)*t0+fxn*t1)
       enddo
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end ratwei

      end subroutine ratwei
