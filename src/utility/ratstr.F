*-----------------------------------------------------------------------
*.specification

      subroutine ratstr (nrid, nz, rtt, nt, rtn, nn, 
     &  rtsa, rtra, rtqa, rtcx, rtrd, rtbr)
      use b2mod_types
      use b2mod_math
      use b2mod_constants
      use b2mod_subsys
      implicit none
#include "emissivity.constants"
*   ..input arguments (unchanged on exit)
      integer nrid, nz, nt, nn
      real (kind=R8) ::
     *   rtt(0:nt), rtn(0:nn)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  rtsa(0:nt,0:nn,0:nz), rtra(0:nt,0:nn,0:nz),
     &  rtqa(0:nt,0:nn,0:nz), rtcx(0:nt,0:nn,0:nz),
     &  rtrd(0:nt,0:nn,0:nz), rtbr(0:nt,0:nn,0:nz)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     RATSTR uses the STRAHL code to prepare tables of rate coefficients
*     for ionisation, recombination, electron energy loss, charge
*     exchange and radiation; as function of temperature and density.
*
*
*  3. description
*
*     This routine employs code from STRAHL /1,2/ to compute a table
*     of atomic rate coefficients.
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
*     The rate coefficients are obtained as in the STRAHL package;
*     a set of routines developed by K. Behringer of Garching and JET
*     (now at the University at Stuttgart).
*     The STRAHL package has been described in Refs. /1,2/; it builds
*     upon Refs. /3-8/ ...
*     For further information, contact K. Behringer at Institut fuer
*     Plasmaforschung, Universitaet Stuttgart, Pfaffenwaldring 31,
*     D-7000 Stuttgart 80, Fed. Rep. Germany.
*
*!!   IMPORTANT NOTE 0: This STRAHL package does not include rate
*!!   coefficients for charge exchange. The value returned by the
*!!   present code is the responsibility of BJB, who has no expertise
*!!   in this matter. The cx rate coefficients supplied here for
*!!   helium might be particularly bad. (See the code for further
*!!   comments.)
*
*!!   IMPORTANT NOTE 1: This STRAHL package does not account for the
*!!   density dependence of the ionisation and energy loss rate
*!!   coefficients. Here, a fudge factor is applied to account for
*!!   the density dependence of the ionisation rate coefficient for
*!!   hydrogen only. The fudge factor is a simple fit to the results
*!!   of Janev in Ref. /9/.
*
*!!   IMPORTANT NOTE 2: This STRAHL package does not account for
*!!   three-body recombination. See the code for further comments and
*!!   for possible treatments.
*
*
*  4. references
*
*     /1/ K. Behringer, Description of the Impurity Transport Code
*     STRAHL, JET-Report R(87)08, 1987.
*     /2/ K. Behringer, ... (journal article?).
*     /3/ W. Lotz, reports IPP 1/62 (1967) and IPP 1/76 (1968).
*     /4/ S. von Goeler et al., Nucl. Fusion 15 (1975) 301.
*     /5/ H. R. Griem, Plasma Spectroscopy, McGraw-Hill (1964).
*     /6/ A. Burgess, Astrophys. J. 141 (1968) 1588.
*     /7/ A. L. Merts et al., report LA-6220-MS (1976).
*     /8/ D. E. Post et al., At. Data Nucl. Data Tables 20 (1977) 397.
*     /9/ R. K. Janev et al., Survey of Atomic Processes in Edge
*     Plasmas, J. Nucl. Mater. 121 (1984) 10-16.
*     /10/ Ya. B. Zeldovich and Yu. P. Raizer: Physics of Shock Waves
*     and High Temperature Hydrodynamic Phenomena, Academic Press, New
*     York, 1966.
*
*
*  5. parameters
*
*     nrid - integer, input.
*     nrid specifies a unit number for formatted input.
*
*     nz - integer, input.
*     nz specifies the nuclear charge of the species for which the rate
*     coefficients are to be evaluated.
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
*     rtsa - (0:nt,0:nn,0:nz) real*8 array, output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtsa(i,j,k) specifies the
*     ionisation rate coefficient at electron temperature rtt(i) and
*     electron density rtn(j) for the process k->k+1. rtsa has
*     dimension m**3/s.
*
*     rtra - (0:nt,0:nn,0:nz) real*8 array, output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtra(i,j,k) specifies the
*     recombination rate coefficient at electron temperature rtt(i) and
*     electron density rtn(j) for the process k->k-1. rtra has
*     dimension m**3/s.
*
*     rtqa - (0:nt,0:nn,0:nz) real*8 array, output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtqa(i,j,k) specifies the
*     electron heat loss rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for processes starting from charge
*     state k. rtqa has dimension eV*m**3/s.
*
*     rtcx - (0:nt,0:nn,0:nz) real*8 array, output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtcx(i,j,k) specifies the charge
*     exchange rate coefficient at neutral hydrogen temperature rtt(i)
*     and electron density rtn(j) for the process k->k-1. rtcx has
*     dimension m**3/s.
*
*     rtrd - (0:nt,0:nn,0:nz) real*8 array, output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtrd(i,j,k) specifies the
*     line radiation rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for processes starting from charge
*     state k. rtrd has dimension eV*m**3/s.
*
*     rtbr - (0:nt,0:nn,0:nz) real*8 array, output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtbr(i,j,k) specifies the
*     bremsstrahlung radiation rate coefficient at electron temperature 
*     rtt(i) and electron density rtn(j) emitted from charge
*     state k. rtbr has dimension eV*m**3/s.
*
*
*  8. further comments
*
*     The present code calculates charge exchange rate coefficients
*     that are manifestly independent of density.  Furthermore, the
*     STRAHL calculation does not consider the density dependence of
*     the ionisation rate coefficients.
*     We apply a 'fudge factor' in order to represent the density
*     dependence of the ionisation rate coefficient for hydrogen only.
*     The fudge factor attempts to reproduce the data of Janev.
*     These issues should be reconsidered.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nzm, i, j, k, klo, khi, kk, nzi
      real (kind=R8) ::
     *  cm3, rydbrg, zlradf, bremsf, zexp1, zexp2, zexp3,
     &  mass, val, zy1, zy2, zy1c, zy2c, zs1, zs2, zs3, zs4, zs5, zs6,
     &  zs7, zz, zzz, cxcoef
      real (kind=R8) ::
     *  t0, t1, t2
      parameter (nzm=99)
      parameter (cm3=1.0e-6_R8, 
     &  zlradf=2.6e-30_R8/ev, bremsf=1.69e-38_R8/ev, 
     1  zexp1 = 1.0_R8/17.0_R8,
     &  zexp2 = 2.0_R8/17.0_R8, zexp3 = 12.0_R8/17.0_R8)
*   ..further input data
      real (kind=R8) ::
     *  d1(0:nzm-1), f1(0:nzm-1), g1(0:nzm-1), d2(0:nzm-1),
     &  f2(0:nzm-1), g2(0:nzm-1), a1(0:nzm-1), b1(0:nzm-1),
     &  c1(0:nzm-1), a2(0:nzm-1), b2(0:nzm-1), c2(0:nzm-1),
     &  x1(0:nzm-1), x2(0:nzm-1), z1(0:nzm-1), z2(0:nzm-1),
     &  qz(0:nzm-1), ze(0:nzm-1)
*   ..procedures
      intrinsic real, sqrt, log, exp, min
      external xertst, xerrab

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('ratstr')
*   ..set cutll, cutlo
      if (.not.b2mod_math_initialised) call b2mod_math_init
*   ..test nz, nt, nn
      call xertst (0.lt.nz.and.nz.le.nzm, 'ratstr--faulty argument nz')
      call xertst (0.le.nt.and.0.le.nn,
     . 'ratstr--faulty argument nt or nn')
*   ..test rtt, rtn
      do i = 0, nt
       call xertst (0.0_R8.lt.rtt(i),
     .  'ratstr--faulty argument array rtt')
      enddo
      do i = 0, nn
       call xertst (0.0_R8.lt.rtn(i),
     .  'ratstr--faulty argument array rtn')
      enddo
*   ..position the input file
      rewind (nrid)
   81 read (nrid,'(//g12.2,i12)',end=82) mass, nzi
      if (nzi.eq.nz) goto 83
      do k = 0, nzi-1
       read (nrid,'(///1x)',end=82)
      enddo
      goto 81
   82 call xerrab ('ratstr--unexpected end of input')
   83 continue
*   ..read data for computation of the rate coefficients
      do k = 0, nz-1
       read (nrid,'(/6g12.4/6g12.4/6g12.4)')
     &   d1(k), f1(k), g1(k), d2(k), f2(k), g2(k),
     &   a1(k), a2(k), b1(k), b2(k), c1(k), c2(k),
     &   x1(k), x2(k), z1(k), z2(k), qz(k), ze(k)
      enddo

* ..compute ionisation rate coefficients rtsa
*   ..note
*     rtsa(i,j,k) will be the rate coefficient at electron temperature
*     rtt(i) and electron density rtn(j), for ionisation from charge
*     state k to charge state k+1. The dimension of rtsa is m**3/s.
*   ..loop over table
      do k = 0, nz-1
       do i = 0, nt
*     ..calculate rtsa(i,0:nn,k)
        zy1 = x1(k)/rtt(i)
        zy2 = x2(k)/rtt(i)
        zy1c = zy1+c1(k)
        zy2c = zy2+c2(k)
        t0 = sei1(zy1)/zy1-b1(k)*sei1(zy1c)/zy1c
        if (a1(k)*z1(k).gt.0.0_R8) then
          t0 = log(a1(k)*z1(k)*t0)-zy1
          if (t0.gt.cutll) then
            zs1 = exp(t0)
          else
            zs1 = 0.0_R8
          endif
        else
          zs1 = 0.0_R8
        end if
        t0 = sei1(zy2)/zy2-b2(k)*sei1(zy2c)/zy2c
        if (a2(k)*z2(k).gt.0.0_R8) then
          t0 = log(a2(k)*z2(k)*t0)-zy2
          zs2 = expu(t0)
        else
          zs2 = 0.0_R8
        end if
        t0 = zs1+zs2
        if (t0.gt.0.0_R8) then
          t0 = log(6.7e-7_R8)+log(t0)-1.5_R8*log(rtt(i))
        else
          t0 = cutll
        end if
        val = expu(t0)
        do j = 0, nn
         rtsa(i,j,k) = expu(log(cm3)+t0)
        enddo
*!!     For the case of nz.eq.1, the values rtsa(,,0) will be corrected
*!!     at the bottom of this routine to account for the dependence on
*!!     density.
       enddo
      enddo
*   ..no ionisation from state nz
      do i = 0, nt
       do j = 0, nn
        rtsa(i,j,nz) = 0.0_R8
       enddo
      enddo

* ..compute recombination rate coefficients rtra
*   ..note
*     rtra(i,j,k) will be the rate coefficient at electron temperature
*     rtt(i) and electron density rtn(j), for recombination from charge
*     state k to charge state k-1. The dimension of rtra is m**3/s.
*   ..no recombination from state 0
      do i = 0, nt
       do j = 0, nn
        rtra(i,j,0) = 0.0_R8
       enddo
      enddo
*   ..contribution from valence shell
      rydbrg=(me/(2.0_R8*qe))*(qe**2/(2.0_R8*const_h*eps0))**2
      do k = 1, nz
       zs2 = k**2*(2*qz(k-1)**2-ze(k-1))/qz(k-1)**3
       do i = 0, nt
        zy1 = x1(k-1)/rtt(i)
        zs1 = sqrt(rydbrg/rtt(i))*zs2*zy1*sei1(zy1)
        do j = 0, nn
         rtra(i,j,k) = cm3*2.6e-14_R8*zs1
        enddo
       enddo
      enddo
*   ..contribution from excited levels
      do k = 1, nz
       do i = 0, nt
        zs2 = 2*k**2*sqrt(rydbrg/rtt(i))
        do j = 0, nn
*        (lower limit: ...)
         klo = nint(qz(k-1))
*        (upper limit: collision limit, stop at klo+20)
         khi = int( min (klo+20.0_R8,
     &     126.0_R8*real(k)**zexp3*(rtt(i)/rydbrg)**zexp1/
     &     (rtn(j)*cm3)**zexp2) )
         do kk = klo, khi
          zy1 = (real(k)/(kk+1))**2*rydbrg/rtt(i)
          rtra(i,j,k) = rtra(i,j,k)+
     &      cm3*2.6e-14_R8*(zs2/(kk+1))*zy1*sei1(zy1)
         enddo
        enddo
       enddo
      enddo
*   ..contribution from di-electronic recombination
*     (exited stages of the recombining ion are important)
      do k = 1, nz-1
       zy1 = d1(k)/(rydbrg*(k+1))
       zy2 = d2(k)/(rydbrg*(k+1))
       if (0.0_R8.le.g1(k).or.nz-2.le.k) then
        zs4 = 1.0_R8+0.105_R8*zy1+0.015_R8*zy1**2
       else
        zs4 = 2.0_R8+0.42_R8*zy1+0.06_R8*zy1**2
       endif
       if (0.0_R8.le.g2(k).or.nz-2.le.k) then
        zs5 = 1.0_R8+0.105_R8*zy2+0.015_R8*zy2**2
       else
        zs5 = 2.0_R8+0.42_R8*zy2+0.06_R8*zy2**2
       endif
       zs1 = f1(k)*sqrt(d1(k))/zs4
       zs2 = f2(k)*sqrt(d2(k))/zs5
       zs3 = 6.5e-10_R8*(k+1)**2*sqrt(k/(k**2+13.4_R8))
       zzz = 1.0_R8+0.015_R8*k**3/(k+1)**2
       do i = 0, nt
        do j = 0, nn
         zz = (1.51e+17_R8/(rtn(j)*cm3)*k**6*sqrt(rtt(i)))**
     1 (1.0_R8/7.0_R8)
         if (0.0_R8.le.g1(k).or.nz-2.le.k) then
          zs6 = zs1/(1.0_R8+200.0_R8/zz)
         else
          zs6 = zs1/(1.0_R8+666.7_R8/((k+1)*zz)**2)
         endif
         if (0.0_R8.le.g2(k).or.nz-2.le.k) then
          zs7 = zs2/(1.0_R8+200.0_R8/zz)
         else
          zs7 = zs2/(1.0_R8+666.7_R8/((k+1)*zz)**2)
         endif
         if (zs7.gt.0.0_R8) then
           t0 = log(zs7)-d2(k)/(rtt(i)*zzz)
           t0 = expu(t0)
         else
           t0 = 0.0_R8
         endif
         if (zs6.gt.0.0_R8) then
           t1 = log(zs6)-d1(k)/(rtt(i)*zzz)
           t1 = expu(t1)
         else
           t1 = 0.0_R8
         end if
         if (t0+t1.gt.0.0_R8) then
           t2 = log(cm3*(zs3/sqrt(rtt(i)**3)))+log(t0+t1)
         else
           t2 = cutll
         end if
         rtra(i,j,k) = rtra(i,j,k)+expu(t2)
        enddo
       enddo
      enddo

* ..compute energy loss rate coefficients rtqa and radiation rate rtrd
*   ..note
*     rtqa(i,j,k) and rtrd(i,j,k) will be the rate coefficients at 
*     electron temperature rtt(i) and electron density rtn(j), for 
*     electron energy loss and radiation from processes on charge 
*     state k, respectively. 
*     The dimensions of rtqa, rtrd and rtbr are eV*m**3/s.
*   ..loop over table
      do k = 0, nz
       do i = 0, nt
        if (k.lt.nz) then
         zs1 = f1(k)*abs(g1(k))*expu(-d1(k)/rtt(i))
         zs2 = f2(k)*abs(g2(k))*expu(-d2(k)/rtt(i))
         if (zs1+zs2.gt.0.0_R8) then
           t0 = log(zlradf)+log(zs1+zs2)-0.5_R8*log(rtt(i))
           t0 = expu(t0)
         else
           t0 = 0.0_R8
         endif
         val = t0 + bremsf*k**2*sqrt(rtt(i))
         do j = 0, nn
          rtqa(i,j,k) = x1(k)*rtsa(i,j,k)+val
          rtrd(i,j,k) = t0
          rtbr(i,j,k) = bremsf*k**2*sqrt(rtt(i))
         enddo
        else
         val = bremsf*k**2*sqrt(rtt(i))
         do j = 0, nn
          rtqa(i,j,k) = val
          rtrd(i,j,k) = 0.0_R8
          rtbr(i,j,k) = val
         enddo
        endif
       enddo
      enddo

* ..compute charge exchange rate coefficients rtcx
*   ..note
*!!   The following code for rtcx is not part of STRAHL.
*     rtcx(i,j,k) will be the rate coefficient for charge exchange
*     recombination from charge state k to charge state k-1 by
*     collision with neutral hydrogen, at normalised temperature
*     rtt(i) and density rtn(j). precisely, rtt(i) is the temperature
*     per atomic mass of the neutral hydrogen, and rtn(j) is the
*     electron density. (rtn is ignored presently.) The dimension of
*     rtcx is m**3/s.
*     A good formula is yet to be supplied. Try perhaps a universal
*     formula of the type rtcx(i,j,k) = c0*rtt(i)**c1*real(k)**c2,
*     with corrections for some reactions.
*     The following assignment is temporary; it comes from a simple
*     fit to results of Puiatti et al. (in: Plasma Physics, Vol. 23
*     (1981) 1075-1092). This fit is the responsibility of BJB, who
*     has no expertise in this matter.
*     The universal formula used here is
*       rtcx = 2.0e-15*E**0.35*Z**1.14
*     (except for hydrogen, see the code).
*     This was obtained from figure 3 of Puiatti, as follows.
*     The figure suggests the approximation
*       sigma = 1.3e-19*E**(-0.15)*Z**1.14
*     give or take 50 percent; preferably restrict to Z.gt.2 and
*     E.lt.2.0e4(ev/amu).
*     If one now sets velocity = sqrt(energy/mass), remembering that
*     E is expressed in eV/amu, then velocity = E**0.5*sqrt(ev/mp).
*     Using ev=1.6e-19 and mp=1.67e-27 then approximately
*       rtcx = 1.3e-15*E**0.35*Z**1.14  (not correct).
*     The actual constant employed is larger, because one computes
*     rtcx = ave(sigma*velocity) by averaging over a Maxwellian
*     distribution. This yields approximately
*       rtcx = 2.0e-15*E**0.35*Z**1.14
*     The error in this formula can be assessed from Puiatti.
*     It is quite bad for Z=2.
*     Note for future reference the work by Janev and Katsonis, in
*     Nuclear Fusion, vol. 27 (1987) pp. 1493ff.
*!!   Reconsider this code, in particular for Z=1 and Z=2.
*   ..no charge exchange recombination from state 0
      do i = 0, nt
       do j = 0, nn
        rtcx(i,j,0) = 0.0_R8
       enddo
      enddo
*   ..simple formula for higher charge states
      do k = 1, nz
*    ..set coefficient cxcoef
*      (cxcoef depends in principle on k and on nz)
       if (nz.eq.1) then
*       (Of course, now also (k.eq.1). This is the case of charge
*       exchange between H0 and H1, a resonant process for which the
*       cross-section is much larger than for H0 upon any other once
*       charged atom).
        cxcoef = 1.3e-14_R8
       else
        cxcoef = 2.0e-15_R8*real(k)**1.14_R8
       endif
       do i = 0, nt
        do j = 0, nn
         rtcx(i,j,k) = cxcoef*rtt(i)**0.35_R8
        enddo
       enddo
      enddo

* ..fudges and fixes by BJB, not part of STRAHL
*   ..density-dependence of rtsa
*!!   We incorporate a fudge factor to take into account the dependence
*!!   on electron density in the ionisation rate coefficient for
*!!   neutral hydrogen only.  Hopefully a more precise treatment will
*!!   become part of STRAHL, in which case the following code should be
*!!   removed.
      if (nz.eq.1) then
       do i = 0, nt
        do j = 0, nn
         rtsa(i,j,0) = rtsa(i,j,0)*hydfud(rtt(i),rtn(j))
        enddo
       enddo
      endif
*     Why is that code placed here at the bottom and not where it
*     belongs, i.e. where rtsa is first computed?  Well, rtsa also
*     occurs as a factor in a term in the energy loss rate coefficient
*     rtqa, but the density dependence of rtqa is in fact much
*     weaker than that of rtsa. We want to modify the STRAHL code only
*     for the most important effect, which is the density dependence
*     of rtsa (and that only for hydrogen). We hope that experts will
*     provide proper code for the other cases.
*   ..three-body processes into rtra
*!!   Estimates such as the Zeldovich formula below show that this
*!!   must be an important process below a few eV at realistic density.
*!!   In order to maintain Saha equilibrium as a valid uniform solution
*!!   it is important to address at the same time multistage ionisation
*!!   processes, i.e. the density dependence of the ionisation rate
*!!   coefficients. It is also important to address the electron energy
*!!   gain associated with the three-body process. At present we do not
*!!   know enough to incorporate all this into the rate coefficients.
*!!   We hope that experts will provide proper code.
*     The following code segment (commented out) follows Zeldovich and
*     the NRL formulary in applying a simple formula that is supposed
*     to be valid for recombination of fully stripped ions. The formula
*     is applied below even in the case that the ion is not fully
*     stripped.  The numerical coefficient may be wrong, for all we
*     know it may be very wrong.
*$$$  do k = 1, nz
*$$$   do i = 0, nt
*$$$    do j = 0, nn
*$$$     rtra(i,j,k) = rtra(i,j,k)+8.75e-39_R8*k**3*rtn(j)/rtt(i)**4.5_R8
*$$$    enddo
*$$$   enddo
*$$$  enddo
*     The following code segment (commented out) shows an alternative
*     treatment of three-body recombination, relying on detailed
*     balance and on Saha equilibrium. We do not know how good or how
*     bad is the approximation to Saha equilibrium that is used below.
*     We also do not know if the following code is better or worse than
*     what is actually done above. In practice the following code
*     segment suffers from the problem that rtsa may have underflowed
*     precisely when the contribution to rtra should be large.
*$$$  do k = 1, nz
*$$$   do i = 0, nt
*$$$    do j = 0, nn
*$$$     rtra(i,j,k) = rtra(i,j,k)+rtsa(i,j,k-1)*rtn(j)*cm3/
*$$$ &     (3.0e21_R8*rtt(i)**1.5_R8*exp(max(-660.0_R8,-x1(k-1)/rtt(i))))
*$$$    enddo
*$$$   enddo
*$$$  enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end ratstr

      contains

      function hydfud(temp,dens)
*!!   (hydfud is the fudge factor for the density dependence in the
*!!   hydrogen ionisation rate coefficient; see the code.)
      implicit none
      real(kind=R8) :: hydfud, temp, dens
      hydfud = exp(4.4_R8*(dens/1.0e20_R8)**0.15_R8*temp**(-0.60_R8))
      return
      end function hydfud

      function sei1(zz)
*     (sei1(z) is an approximation to the scaled exponential integral
*     exp(z)*E1(z))
      implicit none
      real(kind=R8) :: sei1, zz
      sei1 = log(1.0_R8+1.0_R8/zz)-0.4_R8/(1.0_R8+zz)**2
      return
      end function sei1

      end subroutine ratstr
