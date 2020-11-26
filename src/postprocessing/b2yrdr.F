*-----------------------------------------------------------------------
*.specification

      subroutine b2yrdr (ninp, nout)
      use b2mod_types
      use b2mod_version
      use b2mod_b2cmrc
      use b2mod_ppout
      use b2mod_ipmain
      implicit none
#     include "DIMENSIONS.F"
*   ..input arguments (unchanged on exit)
      integer ninp(0:3), nout(0:1)
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YRDR is a postprocessor for B2. It reads and displays a table
*     of atomic rate coefficients.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine displays a table of atomic rate coefficients for the
*     processes of ionisation, recombination, heat loss, and charge
*     exchange. The table is in the form required for use by the B2
*     code /1,2/. It is read from unit ninp(1) and is stored in the
*     common block /b2cmrc/.
*
*     It is immaterial to the present routine how the rate coefficients
*     have been obtained; they may, for instance, have been obtained
*     with use of the ADPAK package or the STRAHL package.
*
*     The ADPAK package was developed by R. A. Hulse of PPL, Princeton,
*     and is described in Refs. /3,4/. For further information, contact
*     R. A. Hulse at Plasma Physics Laboratory, Princeton University,
*     P. O. Box 451, Princeton NJ 08543 (hulse@pppl.gov).
*     The STRAHL package was developed by K. Behringer of IPP Garching
*     and JET (now at the University at Stuttgart). It is described in
*     Refs. /5,6/. For further information, contact K. Behringer at
*     Institut fuer Plasmaforschung, Universitaet Stuttgart,
*     Pfaffenwaldring 31, D-7000 Stuttgart 80, Fed. Rep. Germany
*     (behringe@ipf.uni-stuttgart.de).
*
*     Rate coefficients for hydrogen may instead have been evaluated
*     using tabulated data obtained from D. K. Stotler of PPPL and based
*     on a code developed originally by J. C. Weisheit /7,8/. For
*     further information on these rate coefficients please contact
*     D. K. Stotler at PPPL (postal address above or stotler@pppl.gov).
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
*     ninp - (0:1) integer array, input.
*     ninp specifies the input unit numbers, as follows.
*     ninp(0): formatted; provides output control.
*     ninp(1): un*formatted; provides the table of rate coefficients.
*     ninp(2): un*formatted; provides default physics parameters.
*     ninp(3): un*formatted; provides the geometry.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides printed summary of rate coefs.
*     nout(1): graphical: provides graphical summary of rate coefs.
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer ndismx
      parameter (ndismx=20)
      integer nidis, njdis, idis(0:ndismx-1), jdis(0:ndismx-1)
      integer idum(0:9)
      integer i, is, is1, nx, ny, ns, ns1
      character id*8, lblrc*120, lblcp*120, cnamip*80, cvalip*80
*   ..procedures
      logical streql
      external subini, subend, xertst, xerrab, cfruch
      external cfruin, streql, ipgeti, gacwk, gdawk
      external b2rurc, b2rucp, b2wfr0, b2wfr1, b2wfrc, b2wgrc

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
*     and izhi(inuc) is the highest atomic charge considered in B2 for
*     that nuclear species.
*     rtns equals the variable ns used elsewhere in B2.
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
*     zamin, zamax and zn used elsewhere in B2.
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
*     dimension ev*m**3/s.
*      rtlcx - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlcx(i,j,k) specifies
*     logu(rtrcx), where rtrcx is the charge exchange rate coefficient
*     at neutral hydrogen temperature/amu rtt(i) and electron density
*     rtn(j), for the process k->k-1. rtrcx has dimension m**3/s.
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
*     electron temperature rtt(i), electron density rtn(j), 
*     to reach superstage k. rtlpt has dimension eV. 
*      rtlpi - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlpi(i,j,k) specifies
*     the effective ionisation potential at
*     electron temperature rtt(i), electron density rtn(j), 
*     to reach superstage k from superstage k-1. rtlpi has dimension eV. 
*
*     Description of some local variables:
*
*     lblrc - character*120.
*     lblrc identifies the instance of the program employed to prepare
*     the table of atomic rate coefficients.
*
*     ndismx - integer constant.
*     ndismx is the declared dimension bound for arrays idis and jdis.
*
*     nidis, njdis - integer.
*     nidis and njdis are the effective dimensions of the arrays idis
*     and jdis.
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
      call subini ('b2yrdr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.
     &             1.le.ninp(3).and.1.le.nout(0).and.1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..obtain version numbers
      call cfverr (ninp(1), b2frates_version)
      call cfverr (ninp(2), b2fparam_version)
      call cfverr (ninp(3), b2fgmtry_version)
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(3), 2, idum, 'nx,ny')
      nx = idum(0)
      ny = idum(1)
      call xertst (0.le.nx.and.0.le.ny,
     &  'faulty input nx, ny from geometry file')
*   ..obtain ns from parameters file
      call cfruin (ninp(2), 1, idum, 'ns')
      ns = idum(0)
      call xertst (1.le.ns, 'faulty input ns from parameters file')
*   ..read default physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, nx, ns)
*   ..read label
      call cfruch (ninp(1), 120, lblrc, 'label')
      call get_AM_flag (lblrc, b2frates_flag)
*   ..write label
      write (nout(0),'(/4x,a/4x,a)') lblrc(1:60), lblrc(61:120)
*   ..read output control
      read (ninp(0),'(a8)',err=91) id
      call xertst (streql(id,'*nijdis'),
     &  'faulty input instead of *nijdis, id='//id)
      read (ninp(0),*,err=92) nidis, njdis
      call xertst (0.le.nidis.and.nidis.le.ndismx,
     &  'faulty data for nidis')
      call xertst (0.le.njdis.and.njdis.le.ndismx,
     &  'faulty data for njdis')
      read (ninp(0),'(a8)',err=93) id
      call xertst (streql(id,'*idis'),
     &  'faulty input instead of *idis, id='//id)
      if (1.le.nidis) then
       read (ninp(0),*,err=94) (idis(i),i=0,nidis-1)
      endif
      read (ninp(0),'(a8)',err=95) id
      call xertst (streql(id,'*jdis'),
     &  'faulty input instead of *jdis, id='//id)
      if (1.le.njdis) then
       read (ninp(0),*,err=96) (jdis(i),i=0,njdis-1)
      endif
*   ..read and echo code internal parameters
      write (nout(0),'(/2x,a)') 'non-default internal parameters:'
    1 continue
       read (ninp(0),*,end=2,err=97) cnamip, cvalip
       write (nout(0),'(4x,a,2x,a)') cnamip, trim(cvalip)
       call ipsetc (cnamip, cvalip)
       goto 1
    2 continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'
*   ..check if new bundling or species are required and allocate ns-dependent arrays
      ns1=ns
      call ipgeti ('b2yrdr_ns', ns1)
      call xertst (0.lt.ns1, 'faulty input ns1')
      call xertst (DEF_NSD.ge.ns1,
     . 'increase size of DEF_NSD in DIMENSIONS.F')
*   ..read atomic rate data
      write(*,*) 'ns1 =',ns1
      call alloc_b2mod_ppout(nx,ny,ns1)
      call read_b2mod_ppout('rates', ninp(1), b2frates_flag)
      call b2rurc (ninp(1), b2frates_version)
*   ..check input data
      do i = 0, nidis-1
       call xertst (0.le.idis(i).and.idis(i).lt.rtnt,
     &   'faulty data for idis')
      enddo
      do i = 0, njdis-1
       call xertst (0.le.jdis(i).and.jdis(i).lt.rtnn,
     &   'faulty data for jdis')
      enddo
*   ..write credit
*     (Look for 'adpak', 'strahl' or 'adas' as a substring of the label to
*     decide what credit to write. This is not a rigorous test!)
      select case (b2frates_flag)
      case ('adpak')
       call b2wfr0 (nout(0))
      case ('strahl')
       call b2wfr1 (nout(0))
      case ('adas')
       call b2wfr2 (nout(0))
      case ('amns')
       call b2wfr3 (nout(0))
      case default
       write (nout(0),'(/a/a)')
     &   ' Please check the origin of these rate coefficients;'//
     &   ' this program is unable to determine if they stem from',
     &   ' the ADPAK package, the STRAHL package, ADAS, AMNS,'//
     &   ' or from some other source.'
      end select
*   ..write rtns, rtzmin, rtzmax, rtzn
      write (nout(0),'(/2x,a,i3)')
     &  'number of species in this table: ',rtns
      write (nout(0),'(3x,a,i3,a,f6.2,a,f6.2,a,f6.2)')
     &  ('species ',is,':  zamin =',rtzmin(is),
     &  ' ,  zamax =',rtzmax(is),
     &  ' ,  zn =',rtzn(is), is=0,rtns-1)
*   ..write abscissae
      write (nout(0),'(1p/2x,a/(12e11.3))')
     &  'temperature coordinates: rtt(0:nt)', (rtt(i),i=0,rtnt)
      write (nout(0),'(1p/2x,a/(12e11.3))')
     &  'density coordinates:     rtn(0:nn)', (rtn(i),i=0,rtnn)

* ..display rate coefficients
*   ..activate graphics on nout(1)
      call gacwk (nout(1))
*   ..initialise species index
      is = 0
*   ..loop over nuclear species
   10 if (is.ne.rtns) then
*    ..find range of this nuclear species
       is1 = is+1
   20  if (is1.ne.rtns) then
         if (rtzn(is).eq.rtzn(is1).and.rtzmax(is).lt.rtzmin(is1)) then
         is1 = is1+1
         goto 20
        endif
       endif
*    ..display table
       call b2wfrc (nout(0), is1-is, nidis, idis, njdis, jdis, rtt,
     &   rtnt, rtn, rtnn, rtzmin(is), rtzmax(is), rtzn(is),
     &   fix_recomb, zmax_recomb,
     &   rtlsa(0,0,is), rtlra(0,0,is), rtlqa(0,0,is), rtlcx(0,0,is),
     &   rtlza(0,0,is), rtlz2(0,0,is), rtlpt(0,0,is), rtlpi(0,0,is))
       call b2wgrc (nout(0), is1-is, nidis, idis, njdis, jdis, rtt,
     &   rtnt, rtn, rtnn, rtzn(is), fix_recomb, zmax_recomb,
     &   rtlsa(0,0,is), rtlra(0,0,is), rtlqa(0,0,is), rtlcx(0,0,is),
     &   rtlza(0,0,is), rtlz2(0,0,is), rtlpt(0,0,is), rtlpi(0,0,is))
*    ..increment and loop
       is = is1
       goto 10
      endif
*   ..deactivate graphics on nout(1)
      call gdawk (nout(1))

* ..return
      call dealloc_b2mod_ppout
      call dealloc_b2mod_b2cmrc
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read heading *nijdis')
   92 call xerrab ('error trying to read nidis, njdis data')
   93 call xerrab ('error trying to read heading *idis')
   94 call xerrab ('error trying to read idis data')
   95 call xerrab ('error trying to read heading *jdis')
   96 call xerrab ('error trying to read jdis data')
   97 call xerrab
     1 ('error trying to read non-default internal parameters')

*-----------------------------------------------------------------------
*.end b2yrdr

      end subroutine b2yrdr
