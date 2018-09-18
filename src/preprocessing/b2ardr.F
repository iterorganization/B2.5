*-----------------------------------------------------------------------
*.specification

      subroutine b2ardr (ninp, nout)
      use b2mod_types
      use b2mod_version
      use b2mod_b2cmrc
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:4), nout(0:1)
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2ARDR computes and writes to an external file a table of atomic
*     rate coefficients for b2. The coefficients are obtained using
*     either the ADPAK or the STRAHL set of routines.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine prepares a table of atomic rate coefficients for the
*     processes of ionisation, recombination, heat loss, and charge
*     exchange. The table is in the form required for use by the b2
*     code /1,2/. It is stored in the common block /b2cmrc/ and is
*     written to unit nout(1) for transfer into a b2 calculation.
*
*     The rate coefficients are obtained using either the ADPAK package
*     or the STRAHL package.
*     The ADPAK package was developed by R. A. Hulse of PPL, Princeton,
*     and is described in Refs. /3,4/. For further information, contact
*     R. A. Hulse at Plasma Physics Laboratory, Princeton University,
*     P. O. Box 451, Princeton NJ 08543 (hulse@pppl.gov).
*     The STRAHL package was developed by K. Behringer of IPP Garching
*     and JET (now at the University at Stuttgart).  It is described in
*     Refs. /5,6/. For further information, contact K. Behringer at
*     Institut fuer Plasmaforschung, Universitaet Stuttgart,
*     Pfaffenwaldring 31, D-7000 Stuttgart 80, Fed. Rep. Germany
*     (behringe@ipf.uni-stuttgart.de).  [Since moved back to IPP].
*
*     Selected rate coefficients for hydrogen will be re-evaluated
*     using tabulated data obtained from D. K. Stotler of PPPL and based
*     on a code developed originally by J. C. Weisheit /7,8/. For
*     further information on these rate coefficients please contact
*     D. K. Stotler at PPPL (postal address above or stotler@pppl.gov).
*
*     A smoothing procedure may be applied to the ionisation rate
*     coefficients to soften their steep dependence on temperature at
*     low temperatures. This smoothing is controlled by an input
*     parameter.
*
*
*  4. references
*
*     /1/ B. J. Braams, Computational Studies in Tokamak Equilibrium
*     and Transport, thesis, State University at Utrecht, June 1986.
*     /2/ B. J. Braams, A Multi-Fluid Code for Simulation of the Edge
*     Plasma in Tokamaks, NET Report Nr 68, January 1987.
*     /3/ R. A. Hulse, Numerical Studies of Impurities in Fusion
*     Plasmas, Nuclear Technol./Fusion, vol.3 (1983) pp.259-272.
*     /4/ R. A. Hulse, !! (detailed description?).
*     /5/ K. Behringer, Description of the Impurity Transport Code
*     STRAHL, JET Report R(87)08, 1987.
*     /6/ K. Behringer, !! (journal article?).
*     /7/ Jon C. Weisheit, Recombination in Dense Plasmas, J. Phys. B
*     (Atom. Molec. Phys.) 8 (1975) 2556.
*     /8/ R. K. Janev, D. E. Post, W. D. Langer, K. Evans, D. B.
*     Heifetz and J. C. Weisheit, Survey of Atomic Processes in Edge
*     Plasmas, J. Nuclear Mater. 121 (1984) 10-16.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:2) integer array, input.
*     ninp specifies the input unit numbers, as follows.
*     ninp(0): formatted; provides species parameters.
*     ninp(1): formatted; provides input data for STRAHL calculation.
*     ninp(2): formatted; provides input data from Weisheit code.
*     ninp(3): un*formatted; provides default physics parameters.
*     ninp(4): un*formatted; provides the geometry.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides printed summary of rate coefs.
*     nout(1): un*formatted; provides the table of rate coefficients.
*
*
*  8. further comments
*
*     The routine computes and writes out the logarithm of the atomic
*     rate coefficients in the form rl=logu(r). The function logu is a
*     log protected against underflow: if 0.lt.r then logu(r)=log(r),
*     if r.lt.0 then logu(r) is undefined, and if r.eq.0 then
*     logu(r) is a large negative number so that exp(logu(r)) will
*     surely underflow.
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer nnucm, nzmax, nidis, njdis
      parameter (nnucm=100, nzmax=92, nidis=4, njdis=4)
      integer idis(0:nidis-1), jdis(0:njdis-1)
      integer nnuc, nz(0:nnucm-1), izlo(0:nnucm-1), izhi(0:nnucm-1)
      integer iyear(0:nnucm-1), npart(0:nnucm-1)
      integer izpmin(0:nzmax), izpmax(0:nzmax)
      character flag*8
      integer i, j, k, nx, ny, ns, is, inuc, fix_cx, rtnt_in, rtnn_in,
     *  no_weisheit, fix_recomb, ns_pardf
      integer kstart, kstart1, kend
      integer idum(0:9)
      integer, save :: no_smoothing = 0                        !srv 02.07.12
      character lblrc*120, lblcp*120, lblgm*120, cnamip*80, cvalip*80
      character line*500
      real (kind=R8) ::
     *  lnzero, tailep, tlo, thi, nlo, nhi, cxcoef
      real (kind=R8), allocatable ::
     &  rtsa(:,:,:), rtra(:,:,:),
     &  rtqa(:,:,:), rtcx(:,:,:),
     &  rtrd(:,:,:), rtbr(:,:,:),
     &  rtz2(:,:,:), rtpt(:,:,:),
     &  rtzz(:,:,:), rtpi(:,:,:),
     &  wrk(:,:,:)
      save fix_cx, no_weisheit, rtnt_in, rtnn_in, fix_recomb
      parameter (lnzero=-1.0e30_R8)
      logical unread
*   ..procedures
      intrinsic log, exp
      logical streql, is_comment
      external subini, subend, xertst, xerrab, streql, cfwuch
      external ratstr, ratadp, ratwei, rattai, ratprt, cfruch
      external b2xxid, b2arin, b2wfr0, b2wfr1, b2wurc, b2rucp
      external cfverr, cfruin, is_comment
*   ..initialisation
      data fix_cx/1/, no_weisheit/0/, rtnt_in/40/, rtnn_in/16/,
     & fix_recomb/0/
*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common block has its outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     /b2cmrc/ contains a table of atomic rate coefficients.
*      rtnt - integer constant.
*     rtnt specifies the number of intervals in the discretisation
*     of the temperature range. The table will contain data for
*     (rtnt+1) values of the temperature variable.
*      rtnn - integer constant.
*     rtnn specifies the number of intervals in the discretisation
*     of the density range. The table will contain data for (rtnn+1)
*     values of the density variable.
*      rtns - integer.
*     rtns specifies the total number of atomic species; i.e.
*       (sum inuc : 0.le.inuc.lt.nnuc : izhi(inuc)-izlo(inuc)+1) ,
*     where inuc is a nuclear species index, izlo(inuc) is the lowest
*     and izhi(inuc) is the highest atomic charge considered in b2 for
*     that nuclear species.
*     rtns equals the variable ns used elsewhere in b2.
*      rtdumm - integer.
*     rtdumm is a dummy integer variable, introduced in order to obtain
*     an even-length sequence of integer variables in the common block.
*      rtzmin, rtzmax, rtzn - (0:rtns-1) real*8 array.
*     For 0.le.is.lt.rtns, rtzmin(is) and rtzmax(is) specify the atomic 
*     charge range, and rtzn(is) specifies the nuclear charge
*     of atomic species (is). 
*     The izhi(inuc)-izlo(inuc)+1 charge states for each nuclear
*     species are taken in natural order, starting at izlo(inuc); the
*     highest state of nuclear species (inuc) is followed by the lowest
*     state of (inuc+1). rtzmin, rtzmax and rtzn are equal to the arrays 
*     zamin, zamax and zn used elsewhere in b2.
*      rtt - (0:rtnt) real*8 array.
*     rtt(0:rtnt) specifies the sequence of values of the temperature
*     variable, expressed in ev, for which the rate coefficients are
*     evaluated. rtt(0:rtnt) will be positive and increasing.
*      rtn - (0:rtnn) real*8 array.
*     rtn(0:rtnn) specifies the sequence of values of the density
*     variable, expressed in m**-3, for which the rate coefficients
*     are evaluated. rtn(0:rtnn) will be positive and increasing.
*      rtlt - (0:rtnt) real*8 array.
*     rtlt(0:rtnt) holds log(rtt(0:rtnt)).
*      rtln - (0:rtnn) real*8 array.
*     rtln(0:rtnn) holds log(rtn(0:rtnn)).
*      rtlsa - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlsa(i,j,k) specifies
*     logu(rtrsa), where rtrsa is the ionisation rate coefficient at
*     electron temperature rtt(i), electron density rtn(j), for the
*     process k->k+1. rtrsa has dimension m**3/s. The function logu is
*     the same as (log) except that logu(0.0) is equal to some huge
*     negative number.
*      rtlra - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlra(i,j,k) specifies
*     logu(rtrra), where rtrra is the recombination rate coefficient at
*     electron temperature rtt(i) and electron density rtn(j), for the
*     process k->k-1. rtrra has dimension m**3/s.
*      rtlqa - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlqa(i,j,k) specifies
*     logu(rtrqa), where rtrqa is the electron heat loss rate
*     coefficient at electron temperature rtt(i) and electron density
*     rtn(j), for processes starting from atomic species k. rtrqa has
*     dimension eV*m**3/s.
*      rtlcx - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlcx(i,j,k) specifies
*     logu(rtrcx), where rtrcx is the charge exchange rate coefficient
*     at neutral hydrogen temperature/amu rtt(i) and electron density
*     rtn(j), for the process k->k-1. rtrcx has dimension m**3/s.
*      rtlrd - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlrd(i,j,k) specifies
*     logu(rtrrd), where rtrrd is the line radiation rate
*     coefficient at electron temperature rtt(i) and electron density
*     rtn(j), for processes starting from atomic species k. rtrrd has
*     dimension eV*m**3/s.
*      rtlbr - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlbr(i,j,k) specifies
*     logu(rtrbr), where rtrbr is the bremsstrahlung radiation rate
*     coefficient at electron temperature rtt(i) and electron density
*     rtn(j), for processes starting from atomic species k. rtrbr has
*     dimension eV*m**3/s.
*      rtlza - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlza(i,j,k) specifies
*     the effective charge state at
*     electron temperature rtt(i), electron density rtn(j), for the
*     superstage k. rtlza has dimension e. 
*      rtlz2 - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlz2(i,j,k) specifies
*     the effective squared charge at
*     electron temperature rtt(i), electron density rtn(j), for the
*     superstage k. rtlz2 has dimension e**2. 
*      rtlpt - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlpt(i,j,k) specifies
*     the cumulative ionisation potential at
*     electron temperature rtt(i), electron density rtn(j), to reach
*     superstage k. rtlpt has dimension eV. 
*      rtlpi - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlpi(i,j,k) specifies
*     the effective ionisation potential at
*     electron temperature rtt(i), electron density rtn(j), to reach
*     superstage k from superstage k-1. rtlpi has dimension eV. 
*
*     Description of some local variables:
*
*     nnucm - integer constant.
*     nnucm is the maximum number of nuclear species.
*
*     nzmax - integer constant.
*     nzmax is the maximum allowed value of the nuclear charge.
*
*     nnuc - integer.
*     nnuc specifies the number of nuclear species.
*
*     nz - (0:nnucm-1) integer array.
*     For inuc in (0:nnuc-1), nz(inuc) specifies the nuclear charge of
*     the nuclear species inuc.
*
*     iyear - (0:nnucm-1) integer array.
*     For inuc in (0:nnuc-1), iyear(inuc) specifies the ADAS year of
*     the nuclear species inuc.
*     (if -1 then ADAS year was not the same for all ADAS files)
*
*     rtyr - (0:ntns-1) integer array.
*     Specifies the ADAS year of the atomic species.
*     (if -1 then ADAS year was not the same for all ADAS files)
*
*     npart - (0:nnucm-1) integer array.
*     Specifies the number of partitions of the homonuclear sequence.
*
*     izlo, izhi - (0:nnucm-1) integer array.
*     For inuc in (0:nnuc-1), (izlo(inuc):izhi(inuc)) specifies the
*     range of atomic charge state represented in b2 for the nuclear
*     species inuc. It will hold that
*       0.le.izlo(inuc).le.izhi(inuc)+1 and izhi(inuc).le.nz(inuc).
*
*     flag - character*8.
*     flag specifies which code is to be used to calculate the rate
*     coefficients. Its value may be 'adpak', 'strahl', 'adas' or 'amns'.
*
*     lblrc - character*120.
*     lblrc identifies the instance of the program employed to prepare
*     the table of atomic rate coefficients.
*
*     lnzero - real*8 constant.
*     lnzero is a large negative real*8 number. It must be representable
*     on every computer on which this code might run or might be
*     compiled, but exp(lnzero) must underflow on every such computer.
*
*     tailep - real.
*     tailep specifies a smoothing parameter for the temperature
*     dependence of some rate coefficients. It must be non-negative and
*     small compared to 1; tailep.eq.0.0 implies no smoothing.
*     See the routine rattai for detailed description.
*
*     nidis, njdis - integer constant.
*     nidis and njdis are dimension bounds for arrays idis and jdis.
*
*     idis - (0:nidis-1) integer array.
*     idis is the selection of values of the temperature index at which
*     output will be produced to display the density dependence of the
*     rate coefficients.
*
*     jdis - (0:njdis-1) integer array.
*     jdis is the selection of values of the density index at which
*     output will be produced to display the temperature dependence of
*     the rate coefficients.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ardr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.
     &             1.le.ninp(3).and.1.le.ninp(4).and.1.le.nout(0).and.
     &             1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..obtain version numbers
      call cfverr (ninp(3), b2fparam_version)
      call cfverr (ninp(4), b2fgmtry_version)
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(4), 2, idum, 'nx,ny')
      nx = idum(0)
      ny = idum(1)
      call xertst (0.le.nx.and.0.le.ny,
     &  'faulty input nx, ny from geometry file')
*   ..obtain ns from parameters file
      call cfruin (ninp(3), 1, idum, 'ns')
      ns_pardf = idum(0)
      call xertst (1.le.ns_pardf, 
     &  'faulty input ns from parameters file')
*   ..read header for geometry
      call cfruch (ninp(4), 120, lblgm, 'label')
*   ..obtain range and species data
      call b2arin (ninp(0), nnucm, nzmax, tlo, thi, nlo, nhi,
     &  nnuc, nz, izlo, izhi, flag, tailep)
*   ..read and echo code internal parameters
      write (nout(0),'(/2x,a)') 'non-default internal parameters:'
    1 continue
       read (ninp(0),'(a)',end=2,err=93) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         write (nout(0),'(4x,a,2x,a)') cnamip, trim(cvalip)
         call ipsetc (cnamip, cvalip)
       endif
       goto 1
    2 continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'
*   ..internal parameters
      call ipgeti('b2ardr_fix_cx', fix_cx)
      call ipgeti('b2ardr_fix_recomb', fix_recomb)
      call ipgeti('b2ardr_no_weisheit', no_weisheit)
      call ipgeti('b2ardr_rtnt', rtnt_in)
      call xertst(0.lt.rtnt_in, 'faulty argument rtnt_in')
      call ipgeti('b2ardr_rtnn', rtnn_in)
      call xertst(0.lt.rtnn_in, 'faulty argument rtnn_in')
      call ipgeti('b2ardr_no_smoothing', no_smoothing)              !srv 02.07.12
*   ..obtain label
      call b2xxid ('b2ar', lblrc(1:60))
      if (streql(flag,'adpak')) then
       lblrc(61:120) = 'ADPAK rate coefficients, courtesy of'//
     &   ' R. A. Hulse'
      else if (streql(flag,'strahl')) then
       lblrc(61:120) = 'STRAHL rate coefficients, courtesy of'//
     &   ' K. Behringer'
      else if (streql(flag,'adas')) then
       lblrc(61:120) = 'ADAS rate coefficients, courtesy of'//
     &   ' H. Summers'
#ifdef ITM
      else if (streql(flag,'amns')) then
       lblrc(61:120) = 'AMNS rate coefficients, courtesy of'//
     &   ' EFDA-TF-ITM'
#endif
      else
       call xerrab ('faulty data flag; flag='//flag)
      endif
*   ..echo label
      write (nout(0),'(/2x,a)') 'label for rate coefficients'
      write (nout(0),'(/4x,a/4x,a)') lblrc(1:60), lblrc(61:120)
*   ..write credit
      if (streql(flag,'adpak')) then
       call b2wfr0 (nout(0))
      else if (streql(flag,'strahl')) then
       call b2wfr1 (nout(0))
      else if (streql(flag,'adas')) then
       call b2wfr2 (nout(0))
      endif
*   ..compute and test ns
      ns = 0
      do inuc = 0, nnuc-1
       ns = ns+(1+izhi(inuc)-izlo(inuc))
      enddo
*   ..allocate needed arrays
      call alloc_b2mod_b2cmrc(rtnt_in, rtnn_in, ns)
      allocate(rtsa(0:rtnt_in,0:rtnn_in,0:nzmax))
      allocate(rtra(0:rtnt_in,0:rtnn_in,0:nzmax))
      allocate(rtqa(0:rtnt_in,0:rtnn_in,0:nzmax))
      allocate(rtcx(0:rtnt_in,0:rtnn_in,0:nzmax))
      allocate(rtrd(0:rtnt_in,0:rtnn_in,0:nzmax))
      allocate(rtbr(0:rtnt_in,0:rtnn_in,0:nzmax))
      allocate(rtzz(0:rtnt_in,0:rtnn_in,0:nzmax))
      allocate(rtz2(0:rtnt_in,0:rtnn_in,0:nzmax))
      allocate(rtpi(0:rtnt_in,0:rtnn_in,0:nzmax))
      allocate(rtpt(0:rtnt_in,0:rtnn_in,0:nzmax))
      allocate(wrk(0:rtnt_in,0:rtnn_in,0:nzmax))
*   ..set idis and jdis
*!!   The following assignment is meant to provide output at nice
*     values of the temperature and density coordinate for a specific
*     choice of table parameters. The assignment may be changed as
*     desired.
      do k = 0, nidis-1
       idis(k) = ((k+1)*rtnt)/(nidis+1)
      enddo
      do k = 0, njdis-1
       jdis(k) = ((2*k+1)*rtnn)/(2*njdis)
      enddo

* ..compute abscissae
*   ..set rtns, rtdumm, rtzmin, rtzmax, rtzn
      rtdumm = 0.0_R8
*   ..compute rtlt(0:rtnt) and rtt(0:rtnt)
      do i = 0, rtnt
       rtlt(i) = ((rtnt-i)*log(tlo)+i*log(thi))/rtnt
       rtt(i) = exp(rtlt(i))
      enddo
*   ..compute rtln(0:rtnn) and rtn(0:rtnn)
      do j = 0, rtnn
       rtln(j) = ((rtnn-j)*log(nlo)+j*log(nhi))/rtnn
       rtn(j) = exp(rtln(j))
      enddo
* ..initialise ADAS year
      iyear=0
      npart=0
      rtyr=0
* ..compute and display rate coefficients
*   ..initialise bundled charge model quantities
      rtzz=-1.0_R8
      rtz2=-1.0_R8
      rtpt=-1.0_R8
      rtpi=-1.0_R8
      unread = .true.
*   ..initialise species index
      rtns = 0
      is = 0
*   ..loop over nuclear species
      do inuc = 0, nnuc-1
*    ..compute table
       if (streql(flag,'adpak')) then
        call ratadp (nz(inuc), rtt, rtnt, rtn, rtnn,
     &    rtsa, rtra, rtqa, rtcx, rtrd, rtbr)
        npart(inuc)=nz(inuc)+1
        kstart = izlo(inuc)
        kend = izhi(inuc)
        do k = izlo(inuc), izhi(inuc)
          rtzmin(is+k-izlo(inuc)) = real(k)
          rtzmax(is+k-izlo(inuc)) = real(k)
          rtzn(is+k-izlo(inuc)) = real(nz(inuc))
          rtzz(0:rtnt,0:rtnn,k-izlo(inuc)) = real(k)
          rtz2(0:rtnt,0:rtnn,k-izlo(inuc)) = real(k)**2
        enddo
        rtns = rtns + izhi(inuc) - izlo(inuc) + 1
       else if (streql(flag,'strahl')) then
        call ratstr (ninp(1), nz(inuc), rtt, rtnt, rtn, rtnn,
     &    rtsa, rtra, rtqa, rtcx, rtrd, rtbr)
        npart(inuc)=nz(inuc)+1
        kstart = izlo(inuc)
        kend = izhi(inuc)
        do k = izlo(inuc), izhi(inuc)
          rtzmin(is+k-izlo(inuc)) = real(k)
          rtzmax(is+k-izlo(inuc)) = real(k)
          rtzn(is+k-izlo(inuc)) = real(nz(inuc))
          rtzz(0:rtnt,0:rtnn,k-izlo(inuc)) = real(k)
          rtz2(0:rtnt,0:rtnn,k-izlo(inuc)) = real(k)**2
        enddo
        rtns = rtns + izhi(inuc) - izlo(inuc) + 1
       else if (streql(flag,'adas')) then
        call ratadas (nout(0), nz(inuc), rtt, rtnt, rtn, rtnn,
     &    rtsa, rtra, rtqa, rtcx, rtrd, rtbr, rtzz, rtz2, rtpt, rtpi,
     &    iyear(inuc), fix_recomb,
     &    npart(inuc), izpmin(0:nz(inuc)), izpmax(0:nz(inuc)))
cxpb  We must account for the case when the user does not want all bundles
        kstart = 0
        kend = npart(inuc)-1
        k = kstart
        do while (izpmax(k).lt.izlo(inuc))
          k = k + 1
        end do
        kstart = k
        k = kend
        do while (izpmin(k).gt.izhi(inuc))
          k = k - 1
        end do
        kend = k
        do k = kstart, kend
         rtzmin(is+k-kstart) = izpmin(k)
         rtzmax(is+k-kstart) = izpmax(k)
         rtzn(is+k-kstart) = real(nz(inuc))
        enddo
        rtns = rtns + kend - kstart + 1
#ifdef ITM
       else if (streql(flag,'amns')) then
        call ratamns (nout(0), nz(inuc), rtt, rtnt, rtn, rtnn,
     &    rtsa, rtra, rtqa, rtcx, rtrd, rtbr, rtzz, rtz2, rtpt, rtpi,
     &    iyear(inuc), fix_recomb,
     &    npart(inuc), izpmin(0:nz(inuc)), izpmax(0:nz(inuc)))
cxpb  We must account for the case when the user does not want all bundles
        kstart = 0
        kend = npart(inuc)-1
        k = kstart
        do while (izpmax(k).lt.izlo(inuc))
          k = k + 1
        end do
        kstart = k
        k = kend
        do while (izpmin(k).gt.izhi(inuc))
          k = k - 1
        end do
        kend = k
        do k = kstart, kend
         rtzmin(is+k-kstart) = izpmin(k)
         rtzmax(is+k-kstart) = izpmax(k)
         rtzn(is+k-kstart) = real(nz(inuc))
        enddo
        rtns = rtns + kend - kstart + 1
#endif
       endif
*    ..correct hydrogen data
       if (nz(inuc).eq.1) then
        if (no_weisheit.eq.0) then
         rewind (ninp(2))
         call ratwei (ninp(2), nz(inuc), rtt, rtnt, rtn, rtnn,
     &    rtsa, rtra, rtqa, rtcx, rtrd, rtbr)
        else
         write(*,*)
     1    'No Weisheit correction of the hydrogen rates was performed'
         write(nout(0),*)
     1    'No Weisheit correction of the hydrogen rates was performed'
        endif
       endif
*    ..DPC: fudge for CX
       if ((fix_cx.eq.4).or.(fix_cx.eq.3.and.nz(inuc).eq.1).or.
     1  (fix_cx.eq.2.and.maxval(rtcx(:,:,0:npart(inuc)-1)).le.1e-40_R8)
     .  .or.
     2  (fix_cx.eq.1.and.maxval(rtcx(:,:,0:npart(inuc)-1)).le.1e-40_R8
     3   .and.nz(inuc).eq.1)) then
         write(*,*) 'CX data ''fixed'' for nuclear charge ',nz(inuc)
         write(nout(0),'(/a,i3)')
     1    'CX data ''fixed'' for nuclear charge ',nz(inuc)
         write(*,*) 'Fix based on the fit in ratstr.F by B.J.Braams'
         write(nout(0),'(a)')
     1    'Fix based on the fit in ratstr.F by B.J.Braams'
         if(nz(inuc).eq.1) then
           write(*,*) 'which is probably OK for hydrogenic species'
           write(nout(0),'(a/)')
     1      'which is probably OK for hydrogenic species'
         else
           write(*,*) 'which is problematic for non-hydrogenic species'
           write(nout(0),'(a/)')
     1      'which is problematic for non-hydrogenic species'
         endif
*   ..no charge exchange recombination from neutral state
         if (rtzmax(is).eq.0.0_R8) then
           do i = 0, rtnt
             do j = 0, rtnn
               rtcx(i,j,kstart) = 0.0_R8
             enddo
           enddo
         endif
*   ..simple formula for higher charge states
*   ..kstart1 selects the first charged state of the sequence
         if (rtzmax(is).eq.0.0_R8) then
           kstart1 = kstart + 1
         else
           kstart1 = kstart
         endif
         do k = kstart1, kend
*    ..set coefficient cxcoef
*      (cxcoef depends in principle on k and on nz)
           if (nz(inuc).eq.1) then
*       (Of course, now also (k.eq.1). This is the case of charge
*       exchange between H0 and H1, a resonant process for which the
*       cross-section is much larger than for H0 upon any other once
*       charged atom).
             cxcoef = 1.3e-14_R8
           else
*       We moved the Z**1.14 dependency from cxcoef to rtcx to handle bundles
             cxcoef = 2.0e-15_R8
           endif
           do i = 0, rtnt
             do j = 0, rtnn
               rtcx(i,j,k) = cxcoef*rtt(i)**0.35_R8*rtzz(i,j,k)**1.14_R8
             enddo
           enddo
         enddo
       endif
*    ..tail smoothing for rtsa
       if(no_smoothing.eq.0) then            !srv 02.07.12 {
         call rattai (rtnt, (rtnn+1)*npart(inuc), tailep, rtt, rtsa,
     &    wrk)
         do k = kstart, kend
          do j = 0, rtnn
           do i = 0, rtnt
            rtsa(i,j,k) = wrk(i,j,k)
           enddo
          enddo
         enddo
*    ..tail smoothing for rtqa
         call rattai (rtnt, (rtnn+1)*npart(inuc), tailep, rtt, rtqa,
     &    wrk)
         do k = kstart, kend
          do j = 0, rtnn
           do i = 0, rtnt
            rtqa(i,j,k) = wrk(i,j,k)
           enddo
          enddo
         enddo
*    ..tail smoothing for rtrd
         call rattai (rtnt, (rtnn+1)*npart(inuc), tailep, rtt, rtrd,
     &    wrk)
         do k = kstart, kend
          do j = 0, rtnn
           do i = 0, rtnt
            rtrd(i,j,k) = wrk(i,j,k)
           enddo
          enddo
         enddo
*    ..tail smoothing for rtbr
         call rattai (rtnt, (rtnn+1)*npart(inuc), tailep, rtt, rtbr,
     &    wrk)
         do k = kstart, kend
          do j = 0, rtnn
           do i = 0, rtnt
            rtbr(i,j,k) = wrk(i,j,k)
           enddo
          enddo
         enddo
       endif                                !srv 02.07.12 }
*      (tail is insignificant for rtra and rtcx)
*    ..display table
       call ratprt (nout(0), nz(inuc), npart(inuc), 
     &   nidis, idis, njdis, jdis,
     &   rtt, rtnt, rtn, rtnn, 
     &   rtsa, rtra, rtqa, rtcx, rtrd, rtbr, rtzz, rtz2, rtpt, rtpi)
*    ..enter table into /b2cmrc/
       do k = kstart, kend
        do j = 0, rtnn
         do i = 0, rtnt
          if (k.lt.kend.and.
     &        rtsa(i,j,k).ne.0.0_R8) then
           rtlsa(i,j,is) = log(rtsa(i,j,k))
          else
           rtlsa(i,j,is) = lnzero
          endif
          if (kstart.lt.k.and.rtra(i,j,k).ne.0.0_R8) then
           rtlra(i,j,is) = log(rtra(i,j,k))
          else
           rtlra(i,j,is) = lnzero
          endif
          if (rtqa(i,j,k).ne.0.0_R8) then
           rtlqa(i,j,is) = log(rtqa(i,j,k))
          else
           rtlqa(i,j,is) = lnzero
          endif
          if (rtrd(i,j,k).ne.0.0_R8) then
           rtlrd(i,j,is) = log(rtrd(i,j,k))
          else
           rtlrd(i,j,is) = lnzero
          endif
          if (rtbr(i,j,k).ne.0.0_R8) then
           rtlbr(i,j,is) = log(rtbr(i,j,k))
          else
           rtlbr(i,j,is) = lnzero
          endif
          if (kstart.lt.k.and.rtcx(i,j,k).ne.0.0_R8) then
           rtlcx(i,j,is) = log(rtcx(i,j,k))
          else
           rtlcx(i,j,is) = lnzero
          endif
          if (rtzz(i,j,k).ge.0.0_R8) then
           rtlza(i,j,is) = rtzz(i,j,k)
          else
           rtlza(i,j,is) = real(k)
          endif
          if (rtz2(i,j,k).ge.0.0_R8) then
           rtlz2(i,j,is) = rtz2(i,j,k)
          else
           rtlz2(i,j,is) = real(k**2)
          endif
          if (rtpt(i,j,k).ge.0.0_R8 .and. rtpi(i,j,k).ge.0.0_R8) then
           rtlpt(i,j,is) = rtpt(i,j,k)
           rtlpi(i,j,is) = rtpi(i,j,k)
          else
           if (unread) then
*        ..read default physics common to obtain masses
             call cfruch (ninp(3), 120, lblcp, 'label')
             call b2rucp (ninp(3), b2fparam_version, nx, ns_pardf)
             call initialise_potentials(ns_pardf)
             unread = .false.
           endif
           rtlpt(i,j,is) = (potmin(is)+potmax(is))/2.0_R8
           rtlpi(i,j,is) = (potimin(is)+potimax(is))/2.0_R8
          endif
         enddo
        enddo
*     ..ADAS year
        rtyr(is) = iyear(inuc)
*     ..increment species index
        is = is+1
       enddo
      enddo

* ..save rate table
*   ..version number
      call cfverw (nout(1), newversion)
      call cfwuch (nout(1), 120, lblrc, 'label')
      if (streql(flag,'adas') ) then
#ifdef JET
        call cfwuin (nout(1), rtns, rtyr, 'ppout')
#else
        if (newversion.ge.'01.001.024') then
          call cfwuin (nout(1), rtns, rtyr, 'ppout')
        endif
#endif
      endif
      call b2wurc (nout(1), newversion)

*   ..deallocate arrays
      call dealloc_b2mod_b2cmrc
      deallocate(rtsa)
      deallocate(rtra)
      deallocate(rtqa)
      deallocate(rtcx)
      deallocate(rtrd)
      deallocate(rtbr)
      deallocate(rtzz)
      deallocate(rtz2)
      deallocate(rtpi)
      deallocate(rtpt)
      deallocate(wrk)
* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   93 call xerrab
     1 ('error trying to read non-default internal parameters')

*-----------------------------------------------------------------------
*.end b2ardr

      end subroutine b2ardr
