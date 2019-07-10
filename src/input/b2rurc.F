*-----------------------------------------------------------------------
*.specification

      subroutine b2rurc (nget, version)
      use b2mod_types
      use b2mod_b2cmrc
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nget
      character*10 version
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2RURC reads a table of atomic rate data from an un*formatted
*     file.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine reads from unit (nget) a table of atomic rate data,
*     covering ionisation, recombination, electron heat loss, 
*     charge exchange and radiation with neutral hydrogen. Unit nget
*     must be connected to an un*formatted file.
*     See the code to discern the required format of the input file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nget - integer, input.
*     nget specifies a unit number for un*formatted input.
*
*     (parameters in common:)
*
*     /b2cmrc/ - output.
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
*     electron temperature rtt(i), electron density rtn(j), 
*     to reach superstage k. rtlpt has dimension eV. 
*      rtlpi - (0:rtnt,0:rtnn,0:rtns-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtns-1), rtlpi(i,j,k) specifies
*     the effective ionisation potential at
*     electron temperature rtt(i), electron density rtn(j), 
*     to reach superstage k from superstage k-1. rtlpi has dimension eV. 
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      logical lt0
      integer k, k1, idum(0:9)
      real (kind=R8) ::
     *   cutlo, lnzero, t0
      real (kind=R8), allocatable :: rtza(:), rtza2(:)
      parameter (lnzero=-1.0e30_R8)
*   ..procedures
      intrinsic log
      logical lnext
      real (kind=R8) ::
     *   machsfr, smax
      external subini, subend, xertst, machsfr, smax, cfrure, cfruin,
     1 cfverr
      lnext(k,k1) = rtzmax(k).lt.rtzmin(k1).and.rtzn(k).eq.rtzn(k1)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2rurc')
*   ..test nget
      call xertst (1.le.nget, 'faulty argument nget')
*   ..set cutlo
      cutlo = machsfr()

* ..read rate coefficients
*   ..read dimensions
      call cfruin (nget, 3, idum, 'rtnt,rtnn,rtns')
      call alloc_b2mod_b2cmrc(idum(0), idum(1), idum(2))
*   ..read species data
      if (version.lt.'01.001.024') then 
        allocate(rtza(0:idum(2)-1))   
        allocate(rtza2(0:idum(2)-1))   
        call cfrure (nget, rtns, rtza, 'rtza')
        call cfrure (nget, rtns, rtzn, 'rtzn')
        call cfrure (nget, rtns, rtza2, 'rtza2')
        rtzmin(0:rtns-1) = rtza(0:rtns-1)
        rtzmax(0:rtns-1) = rtza(0:rtns-1)
      else
        call cfrure (nget, rtns, rtzmin, 'rtzmin')
        call cfrure (nget, rtns, rtzmax, 'rtzmax')
        call cfrure (nget, rtns, rtzn, 'rtzn')
      endif
*   ..read abscissae
      call cfrure (nget, rtnt+1, rtt, 'rtt')
      call cfrure (nget, rtnn+1, rtn, 'rtn')
      call cfrure (nget, rtnt+1, rtlt, 'rtlt')
      call cfrure (nget, rtnn+1, rtln, 'rtln')
*   ..read rate coefficients
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtns, rtlsa, 'rtlsa')
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtns, rtlra, 'rtlra')
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtns, rtlqa, 'rtlqa')
      call cfrure (nget, (rtnt+1)*(rtnn+1)*rtns, rtlcx, 'rtlcx')
      if(version.ge.'01.001.016') then  
        call cfrure (nget, (rtnt+1)*(rtnn+1)*rtns, rtlrd, 'rtlrd')
      else
        rtlrd=lnzero
      endif
      if(version.ge.'01.001.019') then  
        call cfrure (nget, (rtnt+1)*(rtnn+1)*rtns, rtlbr, 'rtlbr')
      else
        rtlbr=lnzero
      endif
      if(version.ge.'01.001.024') then  
        call cfrure (nget, (rtnt+1)*(rtnn+1)*rtns, rtlza, 'rtlza')
        call cfrure (nget, (rtnt+1)*(rtnn+1)*rtns, rtlz2, 'rtlz2')
        call cfrure (nget, (rtnt+1)*(rtnn+1)*rtns, rtlpt, 'rtlpt')
        call cfrure (nget, (rtnt+1)*(rtnn+1)*rtns, rtlpi, 'rtlpi')
      else
        call initialise_potentials(rtns)
        do k = 0, rtns-1
          rtlza(0:rtnt,0:rtnn,k)=(rtzmin(k)+rtzmax(k))/2.0_R8
          rtlz2(0:rtnt,0:rtnn,k)=((rtzmin(k)+rtzmax(k))/2.0_R8)**2
          rtlpt(0:rtnt,0:rtnn,k)=(potmin(k)+potmax(k))/2.0_R8
          rtlpi(0:rtnt,0:rtnn,k)=(potimin(k)+potimax(k))/2.0_R8
        enddo
      endif
      if (version.ge.'03.000.007') then
        call cfruin (nget, 1, idum, 'fix_recomb')
        fix_recomb = idum(0)
      else
        fix_recomb = 0 ! corresponding to the default pre-3.0.7 version
      endif 

* ..test rate coefficients
*   ..test ionisation rate coefficients
      do k = 0, rtns-1
        if (k.eq.rtns-1) then
          lt0 = .true.
        else
          lt0 = .not.lnext(k,k+1)
        endif
*      (the meaning of lt0 is 'species k cannot ionise to k+1')
        if (lt0) then
*       (ionisation rate must be 0)
          t0 = smax((rtnt+1)*(rtnn+1),rtlsa(0,0,k),1)
          call xertst (t0.lt.log(cutlo), 'faulty entry in rtlsa table')
        endif
      enddo
*   ..test recombination rate coefficients
      do k = 0, rtns-1
        if (k.eq.0) then
          lt0 = .true.
        else
          lt0 = .not.lnext(k-1,k)
        endif
*      (the meaning of lt0 is 'species k cannot recombine to k-1')
        if (lt0) then
*       (recombination rate must be 0)
          t0 = smax((rtnt+1)*(rtnn+1),rtlra(0,0,k),1)
          call xertst (t0.lt.log(cutlo), 'faulty entry in rtlra table')
        endif
      enddo
*   ..test charge exchange rate coefficients
      do k = 0, rtns-1
        if (k.eq.0) then
          lt0 = .true.
        else
          lt0 = .not.lnext(k-1,k)
        endif
*      (the meaning of lt0 is 'species k cannot recombine to k-1')
        if (lt0) then
*       (charge exchange rate must be 0)
          t0 = smax((rtnt+1)*(rtnn+1),rtlcx(0,0,k),1)
          call xertst (t0.lt.log(cutlo), 'faulty entry in rtlcx table')
        endif
      enddo

      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2rurc
      end subroutine b2rurc

      subroutine get_AM_flag (label, flag)
      implicit none
      character*120 label
      character*8 flag
      integer i, lt(4)
      logical found
      logical streql
      external streql

      flag = 'unknown'
      found = .false.
      lt(1) = index(label,'ADPAK')
      lt(2) = index(label,'STRAHL')
      lt(3) = index(label,'ADAS')
      lt(4) = index(label,'AMNS')
      do i = 1, 4
        if (lt(i).eq.0) cycle
        if (found) then
          flag = 'multiple'
        else
          select case (i)
          case (1)
            flag = 'adpak'
          case (2)
            flag = 'strahl'
          case (3)
            flag = 'adas'
          case (4)
            flag = 'amns'
          end select
        end if
        found = .true.
      end do
      call xertst(.not.streql(flag,'unknown'),
     . 'Could not identify A&M data package !')
      call xertst(.not.streql(flag,'multiple'),
     . 'Ambiguous multiple A&M data package setting !')

      return
      end
