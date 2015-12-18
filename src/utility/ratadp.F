*-----------------------------------------------------------------------
*.specification

      subroutine ratadp (nz, rtt, nt, rtn, nn, 
     1  rtsa, rtra, rtqa, rtcx, rtrd, rtbr)
      use b2mod_types
      use b2mod_constants
      implicit none
*   ..input arguments (unchanged on exit)
      integer nz, nt, nn
      real (kind=R8) ::
     *  rtt(0:nt), rtn(0:nn)
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
*     RATADP uses the ADPAK code to prepare tables of rate coefficients
*     for ionisation, recombination, electron energy loss, charge
*     exchange, and radiation rate; as function of temperature and density.
*
*
*  3. description
*
*     This routine employs ADPAK /1,2/ to compute a table of atomic
*     rate coefficients.
*
*     rtt specifies temperature or temperature per amu; units are ev.
*     rtn specifies density; units are m**(-3).
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
*     The rate coefficients are obtained with use of the ADPAK package;
*     a set of routines developed by R. A. Hulse of PPL, Princeton. The
*     ADPAK package has been described in Refs. /1,2/; it builds upon
*     work by Post et al. /3/, Hulse et al. /4/, and Bell et al. /5/.
*     For further information, contact R. A. Hulse at Plasma Physics
*     Laboratory, Princeton University, P. O. Box 451, Princeton
*     NJ 08543.
*
*!!   IMPORTANT NOTE: The ADPAK package does not account for the
*!!   density dependence of the ionisation and energy loss rate
*!!   coefficients. Here, a fudge factor is applied to account for
*!!   the density dependence of the ionisation rate coefficient for
*!!   hydrogen only. The fudge factor is a simple fit to the results
*!!   of Janev in Ref. /6/.
*
*
*  4. references
*
*     /1/ R. A. Hulse, "Numerical Studies of Impurities in Fusion
*     Plasmas", Nuclear Technol./Fusion, vol.3 (1983) pp.259-272.
*     /2/ R. A. Hulse, ... (detailed description?).
*     /3/ D. E. Post, R. V. Jensen, C. B. Tarter, W. H. Grasberger,
*     and W. A. Lokke, "Steady-State Radiative Cooling Rates for Low-
*     Density, High-Temperature Plasma", At. Data Nucl. Data Tables
*     vol.20 (1977) pp.397-.
*     /4/ R. A. Hulse, D. E. Post, and D. R. Mikkelsen, "Charge
*     Exchange as a Recombination Mechanism in High-Temperature
*     Plasmas", J. Phys. B, vol.13 (1980) pp.3895-.
*     /5/ K. L. Bell, H. B. Gilbody, J. G. Hughes, A. E. Kingston,
*     and F. J. Smith, "Recommended Data on the Electron Impact
*     Ionization of Light Atoms and Ions", J. Phys. Chem. Ref. Data,
*     vol.12 (1983) pp.891-916
*     /6/ R. K. Janev et al., "Survey of Atomic Processes in Edge
*     Plasmas", J. Nucl. Mater. 121 (1984) 10-16.
*
*
*  5. parameters
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
*     The present code calculates charge-exchange rate coefficients
*     that are manifestly independent of density.  Furthermore, the
*     inner workings of the ADPAK routines do not consider the
*     density dependence of the ionisation rate coefficients.
*     We apply a 'fudge factor' in order to represent the density
*     dependence of the ionisation rate coefficient for hydrogen only.
*     The fudge factor attempts to reproduce the data of Janev.
*     These issues should be re-considered.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nzm, i, j, k
      real (kind=R8) ::
     *  cm3, arbn, rad, temp, dens
      parameter (nzm=99, cm3=1.0e-6_R8, arbn=1.0e18_R8)
*     cm3 is one cm**3, expressed in m**3.
*     arbn is an arbitrary density, expressed in 1/m**3.
*   ..common blocks for adpak
      integer nucz, nspc, nvalnc(100)
      real (kind=R8) ::
     *  apn(100,10), aqn(100,10), sigma(10,10), qn(100,10),
     &  en(100,10), ein(100,5), enm(100,5,10),
     &  fnm(100,5,10), fnn(100), enn(100)
      common /adsdat/ nucz, nspc, apn, aqn, nvalnc, sigma, qn,
     &  en, ein, enm, fnm, fnn, enn
      real (kind=R8) ::
     *  rclion(100), rrarec(100), rdirec(100), cizlos(100),
     &  radrrc(100), raddrc(100), radclx(100), radbrm(100)
      common /adelec/ rclion, rrarec, rdirec, cizlos,
     &  radrrc, raddrc, radclx, radbrm
      real (kind=R8) ::
     *  rcxrec(100), radbcx(100)
      common /adneut/ rcxrec, radbcx
      integer laden, ladtip, ltipok, leci, leciok
      common /adparm/ laden, ladtip, ltipok, leci, leciok
      integer ldrmlt
      real (kind=R8) ::
     *  cdnn(100), cdnm(100),
     &  yhnna, yhnnb, yhnnc, yhnma, yhnmb, yhnmc
      common /addiel/ cdnn, cdnm,
     &  yhnna, yhnnb, yhnnc, yhnma, yhnmb, yhnmc, ldrmlt
      real (kind=R8) ::
     *  cefrac(100), ceavgz, cerad
      common /adceq/ cefrac, ceavgz, cerad
*   ..further interface to adpak
      integer ncxopt, ncxerr, ivunit
      real (kind=R8) ::
     *   anneut(0:0), vneut(0:0)
*   ..procedures
      real (kind=R8) ::
     *   hydfud
      external adset, aderc, adbcxr, subini, subend, xertst
      hydfud(temp,dens) = exp(4.4_R8*(dens/1.0e20_R8)**0.15_R8*
     1 temp**(-0.60_R8))
*!!   (hydfud is the fudge factor for the density dependence in the
*!!   hydrogen ionisation rate coefficient; see the code.)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('ratadp')
*   ..test nz, nt, nn
      call xertst (0.lt.nz.and.nz.le.nzm, 'ratadp--faulty argument nz')
      call xertst (0.le.nt.and.0.le.nn, 
     . 'ratadp--faulty argument nt or nn')
*   ..test rtt, rtn
      do i = 0, nt
       call xertst (0.0_R8.lt.rtt(i),
     .  'ratadp--faulty argument array rtt')
      enddo
      do i = 0, nn
       call xertst (0.0_R8.lt.rtn(i),
     . 'ratadp--faulty argument array rtn')
      enddo
*   ..set some adpak parameters
*     (contact R. A. Hulse for advice)
      laden  = 1
      ladtip = 1
      leci   = 2
      ldrmlt = 1
      ncxopt = 2
      ivunit = 2

* ..compute rate coefficients
*   ..initialisation
      call adset (nz)
*   ..loop over table
      do i = 0, nt
       do j = 0, nn
*     ..compute rates for electron processes
*       (aderc takes arguments: te/kev, ne*cm3)
        call aderc (rtt(i)/1.0e3_R8, rtn(j)*cm3)
*     ..compute rate for charge exchange
*       (adbcxr takes arguments: t0/kev, n0*cm3)
        vneut(0) = rtt(i)/1.0e3_R8
        anneut(0) = arbn*cm3
        call adbcxr (anneut, vneut, 1, ncxopt, ncxerr, ivunit)
*     ..test adpak results
        call xertst (ncxerr.eq.0, 'ratadp--error, ncxerr not 0')
        call xertst (rclion(nz+1).eq.0.0_R8,
     .   'ratadp--error, ionisation from z=nz not 0')
        call xertst (rrarec(1).eq.0.0_R8.and.rdirec(1).eq.0.0_R8,
     .   'ratadp--error, recombination from z=0 not 0')
        call xertst (rcxrec(1).eq.0.0_R8,
     .   'ratadp--error, cx-recombin from z=0 not 0')
*     ..compute rate coefficients
        do k = 0, nz
*        (ionisation rate coefficient; units are m**3/s)
         rtsa(i,j,k) = rclion(k+1)/rtn(j)
*        (recombination rate coefficient; units are m**3/s)
         rtra(i,j,k) = (rrarec(k+1)+rdirec(k+1))/rtn(j)
*        (electron energy loss rate coefficient; eV*m**3/s)
         rad = radclx(k+1)+cizlos(k+1)+radbrm(k+1)+
     &     rtt(i)*ev*rrarec(k+1)+rtt(i)*ev*rdirec(k+1)
         rtqa(i,j,k) = rad/(ev*rtn(j))
*        (charge exchange rate coefficient; units are m**3/s)
         rtcx(i,j,k) = rcxrec(k+1)/arbn
*        (line radiation rate coefficient; eV*m**3/s)
         rad = radclx(k+1)+radrrc(k+1)+raddrc(k+1)
         rtrd(i,j,k) = rad/(ev*rtn(j))
*        (bremsstrahlung radiation rate coefficient; eV*m**3/s)
         rtbr(i,j,k) = radbrm(k+1)/(ev*rtn(j))
        enddo
       enddo
      enddo

* ..density-dependence of rtsa
*!!   (We incorporate a fudge factor to take into account the dependence
*!!   on electron density in the ionisation rate coefficient for
*!!   neutral hydrogen only.  Hopefully a more precise treatment will
*!!   become part of ADPAK, in which case the following code should be
*!!   removed.)
      if (nz.eq.1) then
       do i = 0, nt
        do j = 0, nn
         rtsa(i,j,0) = rtsa(i,j,0)*hydfud(rtt(i),rtn(j))
        enddo
       enddo
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end ratadp

      end subroutine ratadp
