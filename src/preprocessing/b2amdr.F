*-----------------------------------------------------------------------
*.specification

      subroutine b2amdr (ninp, nout)
      use b2mod_types
      use b2mod_version
      use b2mod_b2cmrc_mol
      use b2mod_b2cmpa
      use b2mod_ipmain
      use b2mod_subsys
      use b2mod_math
      use b2mod_constants
#if ( defined(IMAS) && defined(AMNS) )
      use amns_types       ! IGNORE
      use amns_module      ! IGNORE
#elif ( defined(ITM) && defined(GOT_AMNSPROTO) )
      use amns_types       ! IGNORE
      use amns_module      ! IGNORE
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp, nout(0:1)
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AMDR computes and writes to an external file a table of molecular
*     rate coefficients for B2.5. 
*
*
*  3. description (see also routine b2cdca)
*
*     This routine prepares a table of molecular rate coefficients.
*
*
*
*  4. references
*
*     !! TODO
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - integer , input.
*     ninp specifies the input unit numbers, as follows.
*     ninp: formatted; provides species parameters.

*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides printed summary of rate coefs.
*     nout(1): un*formatted; provides the table of rate coefficients.
*
*
*  8. further comments
*
*     The routine computes and writes out the logarithm of the molecular
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
      integer nms_max, nms, ims, N_int
      parameter (nms_max=100)
      integer i, j, k, l, ns, is, js,
     *  rtnt_in, rtnn_in, rtntr_in, rtn_mcoll_in
      integer kstart, kstart1, kend
      integer, allocatable :: is1(:)
      integer idum(0:9)
      integer, save :: no_smoothing = 0              
      integer, save :: ismain0 = 0
      character lblrc*120, lblcp*120, cnamip*80, cvalip*80
      character line*500, ligne*500
      real (kind=R8) ::
     *  lnzero, tailep, tlo, thi, nlo, nhi, cxcoef, am0, am1, am2,
     *  tr_min, tr_max, mcoll_min, mcoll_max
      real (kind=R8), allocatable ::
     &  rt_d(:,:,:),   rt_nd(:,:,:),
     &  rt_di1(:,:,:), rt_ec_m(:,:,:,:),
     &  rt_ec_E(:,:,:),
     &  rt_ic(:,:,:),  rt_dr(:,:,:),
     &  rt_de(:,:,:),  rt_di2(:,:,:)
#if ( defined(IMAS) && defined(AMNS) )
      type (amns_handle_type), save :: amns
#else
      character*5, save :: git_version_amns = '0.0.0'
#endif
      save rtnt_in, rtnn_in
      parameter (lnzero=-1.0e30_R8)
*   ..procedures
      intrinsic log, exp
      logical streql, is_comment
      external xertst, xerrab, streql, cfwuch
      external ratamjhyd_mol
      external ratstr, cfruch
      external b2xxid, b2amin, b2wfr0, b2wfr1, b2wurc, b2rucp
      external cfruin, is_comment, get_ADAS_hash
*   ..initialisation
      data rtnt_in/40/, rtnn_in/32/
*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common block has its outermost declaration in
*     this routine; they need not be preserved between calls.
*
*      rtnt - integer constant.
*     rtnt specifies the number of intervals in the discretisation
*     of the temperature range. The table will contain data for
*     (rtnt+1) values of the temperature variable.
*
*      rtnn - integer constant.
*     rtnn specifies the number of intervals in the discretisation
*     of the density range. The table will contain data for (rtnn+1)
*     values of the density variable.
*
*      rtms - integer.
*     rtms specifies the total number of molecular species.
*
*      rtzn1, rtzn2 - (0:rtnms-1) real*8 array.
*     For 0.le.ims.lt.rtnms, rtzn1(ims) specifies the nuclear charge
*     of the first atom of the molecule (ims). rtzn2 contains the second. 
*
*      rtt - (0:rtnt) real*8 array.
*     rtt(0:rtnt) specifies the sequence of values of the temperature
*     variable, expressed in ev, for which the rate coefficients are
*     evaluated. rtt(0:rtnt) will be positive and increasing.
*
*      rtn - (0:rtnn) real*8 array.
*     rtn(0:rtnn) specifies the sequence of values of the density
*     variable, expressed in m**-3, for which the rate coefficients
*     are evaluated. rtn(0:rtnn) will be positive and increasing.
*
*      rtlt - (0:rtnt) real*8 array.
*     rtlt(0:rtnt) holds log(rtt(0:rtnt)).
*
*      rtln - (0:rtnn) real*8 array.
*     rtln(0:rtnn) holds log(rtn(0:rtnn)).
*
*      rtl_d - (0:rtnt,0:rtnn,0:rtnms-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtnms-1), rtl_d(i,j,k) specifies
*     logu(rt_d), where rt_d is the dissociation rate coefficient at
*     electron temperature rtt(i), electron density rtn(j).
*     rt_d has dimension m**3/s. The function logu is
*     the same as (log) except that logu(0.0) is equal to some huge
*     negative number.
*
*      rtl_nd - (0:rtnt,0:rtnn,0:rtnms-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtnms-1), rtl_nd(i,j,k) specifies
*     logu(rt_d), where rt_d is the rate coefficient for non-dissociative
*     ionization at
*     electron temperature rtt(i), electron density rtn(j).
*     rt_d has dimension m**3/s.
*
*      rtl_di1 - (0:rtnt,0:rtnn,0:rtnms-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtnms-1), rtl_di1(i,j,k) specifies
*     logu(rt_di1), where rt_di1 is the rate coefficient for dissociative
*     ionization (from a molecule) at
*     electron temperature rtt(i), electron density rtn(j).
*     rt_di1 has dimension m**3/s.
*
*      rtl_ec_m - (0:rtnt,0:rtntr,0:rtn_mcoll,0:rtnms-1) real*8 array.
*     For (i,j,l,k) in (0:rtnt,0:rtnn,0:rtn_mcoll,0:rtnms-1), rtl_ec_m(i,j,l,k) specifies
*     logu(rt_ec_m), where rt_ec_m is the momentum rate coefficient for elastic collisions
*     between ions and molecules at 
*     rt_ec_m has dimension m**3/s.
*
*      rtl_ec_E - (0:rtnt,0:rtntr,0:rtnms-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtnms-1), rtl_ec_E(i,j,k) specifies
*     logu(rt_ec_E), where rt_ec_E is the energy rate coefficient for elastic collisions
*     between ions and molecules at 
*     rt_ec_E has dimension m**3/s.
*
*      rtl_ic - (0:rtnt,0:rtnn,0:rtnms-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtnms-1), rtl_ic(i,j,k) specifies
*     logu(rt_ic), where rt_ic is the rate coefficient for ion conversion at
*     ion temperature rtt(i). 
*     The dependency on the electron density is kept in the table structurally, but has no effect for ic.
*     rt_ic has dimension m**3/s.
*
*      rtl_dr - (0:rtnt,0:rtnn,0:rtnms-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtnms-1), rtl_dr(i,j,k) specifies
*     logu(rt_dr), where rt_dr is the rate coefficient for dissociative recombination
*     from a molecular ion at
*     electron temperature rtt(i), electron density rtn(j).
*     rt_dr has dimension m**3/s.
*
*      rtl_de - (0:rtnt,0:rtnn,0:rtnms-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtnms-1), rtl_de(i,j,k) specifies
*     logu(rt_de), where rt_de is the rate coefficient for dissociative excitation
*     from a molecular ion at
*     electron temperature rtt(i), electron density rtn(j).
*     rt_de has dimension m**3/s.
*
*      rtl_di2 - (0:rtnt,0:rtnn,0:rtnms-1) real*8 array.
*     For (i,j,k) in (0:rtnt,0:rtnn,0:rtnms-1), rtl_di2(i,j,k) specifies
*     logu(rt_di2), where rt_di2 is the rate coefficient for dissociative ionization
*     from a molecular ion at
*     electron temperature rtt(i), electron density rtn(j).
*     rt_di2 has dimension m**3/s.
*
*     Description of some local variables:
*
*     nms_max - integer constant.
*     nms_max is the maximum number of nuclear species.
*
*     nms - integer.
*     nms specifies the number of nuclear species.
*
*     lblrc - character*120. !! ????????????
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
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2amdr')
*   ..test ninp, nout
      call xertst (1.le.ninp.and.1.le.nout(0).and.
     &             1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..obtain nms from parameters file
cc      call cfruin (ninp, 1, idum, 'nms')
*   ..obtain range and species data
      call b2amin (ninp, nms_max, tlo, thi, nlo, nhi,
     &  nms, tailep)
*   ..read and echo code internal parameters
      write (nout(0),'(/2x,a)') 'non-default internal parameters:'
    1 continue
       read (ninp,'(a)',end=2,err=93) line
       if (.not.is_comment(line)) then
         ligne = line
         call strip_spaces(ligne)
         if (ligne(1:1).eq.'''') then
           read (line,*) cnamip, cvalip
           write (nout(0),'(4x,a,2x,a)') cnamip, trim(cvalip)
           call ipsetc (cnamip, cvalip)
         end if
       endif
       goto 1
    2 continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'
*   ..internal parameters
      call ipgeti('b2amdr_rtnt', rtnt_in)
      call xertst(0.lt.rtnt_in, 'faulty argument rtnt_in')
      call ipgeti('b2amdr_rtnn', rtnn_in)
      call xertst(0.lt.rtnn_in, 'faulty argument rtnn_in')
      call ipgeti('b2amdr_no_smoothing', no_smoothing)        
      call ipgeti('b2amdr_rtntr', rtntr_in)
      call xertst(0.lt.rtnt_in, 'faulty argument rtntr_in')
      call ipgeti('b2amdr_rtn_mcoll', rtn_mcoll_in)
      call xertst(0.lt.rtn_mcoll_in, 'faulty argument rtn_mcoll_in') 

      call ipgetr('b2amdr_tr_min', tr_min)
      call ipgetr('b2amdr_tr_max', tr_max)
      call ipgetr('b2amdr_mcoll_min', mcoll_min)
      call ipgetr('b2amdr_mcoll_max', mcoll_max)

*   ..obtain label
      call b2xxid ('b2am', lblrc(1:60)) !! TODO

*   ..echo label
      write (nout(0),'(/2x,a)') 'label for rate coefficients'
      write (nout(0),'(/4x,a/4x,a)') lblrc(1:60), lblrc(61:120)

*   ..compute and test nms
      !! TODO

*   ..allocate needed arrays
      call alloc_b2mod_b2cmrc_mol(rtnt_in, rtnn_in, nms,
     .  rtntr_in, rtn_mcoll_in)
      allocate(rt_d(0:rtnt_in,0:rtnn_in,0:nms))
      allocate(rt_nd(0:rtnt_in,0:rtnn_in,0:nms))
      allocate(rt_di1(0:rtnt_in,0:rtnn_in,0:nms))
      allocate(rt_ec_m(0:rtnt_in,0:rtntr_in,0:rtn_mcoll_in,
     .         0:nms))
      allocate(rt_ec_E(0:rtnt_in,0:rtntr_in,0:nms))
      allocate(rt_ic(0:rtnt_in,0:rtnn_in,0:nms))
      allocate(rt_dr(0:rtnt_in,0:rtnn_in,0:nms))
      allocate(rt_de(0:rtnt_in,0:rtnn_in,0:nms))
      allocate(rt_di2(0:rtnt_in,0:rtnn_in,0:nms))

* ..compute abscissae
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
*   ..compute rtltr(0:rtntr) and rttr(0:rtntr). This one on linear scale instead of log
      do j = 0, rtntr
       rttr(j) = ((rtntr-j)*tr_min+j*tr_max)/rtntr
       rtltr(j) = log(rttr(j))
      enddo
*   ..compute rtl_mcoll(0:rtn_mcoll) and rt_mcoll(0:rtn_mcoll). This one on linear scale instead of log
      do j = 0, rtn_mcoll
       rt_mcoll(j) = ((rtn_mcoll-j)*mcoll_min+j*
     .                 mcoll_max)/rtn_mcoll
       rtl_mcoll(j) = log(rt_mcoll(j))
      enddo

* ..compute rate coefficients

*   ..initialise species index
      rtnms = nms
      ims = 0


*   ..loop over nuclear species
      do ims = 0, nms-1
       am1 = 2.0_R8 !! TODO
       am2 = 2.0_R8
*    ..compute table
       write(*,*)
     1    'Molecular rates were taken from Amjuel-Hydhel'
       write(nout(0),*)
     1    'Molecular rates were taken from Amjuel-Hydhel'
       rewind (ninp)
       call ratamjhyd_mol(ninp, nms, am1, am2, rtt, rtnt, rtn, rtnn,
     &  rttr, rtntr, rt_mcoll, rtn_mcoll, 
     &  rt_nd, rt_d, rt_di1, rt_ec_m, rt_ec_E, rt_ic, rt_dr, 
     &  rt_de, rt_di2)
      enddo !! end loop over molecular species

*    ..enter table into /b2cmrc/
      do ims = 0, nms-1
        do j = 0, rtnn
          do i = 0, rtnt
            if (rt_nd(i,j,ims).ne.0.0_R8) then
              rtl_nd(i,j,ims) = log(rt_nd(i,j,ims))
            else
              rtl_nd(i,j,ims) = lnzero
            endif

            if (rt_d(i,j,ims).ne.0.0_R8) then
              rtl_d(i,j,ims) = log(rt_d(i,j,ims))
            else
              rtl_d(i,j,ims) = lnzero
            endif

            if (rt_di1(i,j,ims).ne.0.0_R8) then
              rtl_di1(i,j,ims) = log(rt_di1(i,j,ims))
            else
              rtl_di1(i,j,ims) = lnzero
            endif

            if (rt_ic(i,j,ims).ne.0.0_R8) then
              rtl_ic(i,j,ims) = log(rt_ic(i,j,ims))
            else
              rtl_ic(i,j,ims) = lnzero
            endif

            if (rt_dr(i,j,ims).ne.0.0_R8) then
              rtl_dr(i,j,ims) = log(rt_dr(i,j,ims))
            else
              rtl_dr(i,j,ims) = lnzero
            endif

            if (rt_de(i,j,ims).ne.0.0_R8) then
              rtl_de(i,j,ims) = log(rt_de(i,j,ims))
            else
              rtl_de(i,j,ims) = lnzero
            endif

            if (rt_di2(i,j,ims).ne.0.0_R8) then
              rtl_di2(i,j,ims) = log(rt_di2(i,j,ims))
            else
              rtl_di2(i,j,ims) = lnzero
            endif
          enddo
        enddo
      enddo

cc elastic collisions are treated differently
      do ims = 0, nms-1
        do i = 0, rtnt
          do k = 0, rtntr_in
            if (rt_ec_E(i,k,ims).ne.0.0_R8) then
                rtl_ec_E(i,k,ims) = log(rt_ec_E(i,k,ims))
            else
                rtl_ec_E(i,k,ims) = lnzero
            endif
            do l = 0, rtn_mcoll_in
              if (rt_ec_m(i,k,l,ims).ne.0.0_R8) then
                rtl_ec_m(i,k,l,ims) = log(rt_ec_m(i,k,l,ims))
              else
                rtl_ec_m(i,k,l,ims) = lnzero
              endif
            enddo
          enddo
        enddo
      enddo

* ..save rate table
*   ..version number
cc      call cfverw (nout(1), newversion) !! TODO revive?
cc      call cfwuch (nout(1), 120, lblrc, 'label')
      call b2wurc_mol (nout(1), newversion)

*   ..deallocate arrays
      call dealloc_b2mod_b2cmrc_mol
      deallocate(rt_nd)
      deallocate(rt_d)
      deallocate(rt_di1)
      deallocate(rt_ec_m)
      deallocate(rt_ec_E)
      deallocate(rt_ic)
      deallocate(rt_dr)
      deallocate(rt_de)
      deallocate(rt_di2)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   93 call xerrab
     1 ('error trying to read non-default internal parameters')

*-----------------------------------------------------------------------
*.end b2amdr

      end subroutine b2amdr
