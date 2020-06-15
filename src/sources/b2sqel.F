*-----------------------------------------------------------------------
*.specification

      subroutine b2sqel (nCv, ns, ismain, ev, te, rt, rtw)
      use b2mod_types
      use b2mod_b2cmpa
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns, ismain
      real (kind=R8) ::
     *  ev, te(1:nCv)
      type (B2Rates), intent (inout) :: rt
*   ..output arguments (unspecified on entry)
      type (B2RatesWork), intent (inout) :: rtw

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SQEL computes on all cells of the mesh the rate coefficients
*     for ionisation, recombination, electron heat loss and radiation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine employs the log-log linear representation provided
*     in the arrays rlsa, rlra, rlqa, rlrd and rlbr to evaluate on all cells 
*     of the mesh the rate coefficients for ionisation, recombination, 
*     electron heat loss, and radiation.
*
*
*  8. further comments
*
*     In case a rate coefficient vanishes, the representation of the
*     log of the rate coefficient will provide a large negative value,
*     for which the exponential function should underflow. The present
*     code tests for this condition and explicitly sets those rate
*     coefficients to 0.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv, is
      real (kind=R8) ::
     *   phm0, phm1, phm2, art_rad, tnorm
      real (kind=R8) ::
     *  wrk0(1:nCv), expu
      save ncall, phm0, phm1, phm2, art_rad
      integer, save :: use_eirene=0
*   ..procedures
      intrinsic log, exp, min, max, dim
      external subini, subend, xertst, ipgetr, expu
*   ..initialisation
      data ncall/0/, phm0/1.0_R8/, phm1/1.0_R8/, phm2/1.0_R8/,
     1 art_rad/0.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sqel')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2sqel_phm0', phm0)
       call ipgetr ('b2sqel_phm1', phm1)
       call ipgetr ('b2sqel_phm2', phm2)
       call xertst (0.0_R8.le.min(phm0,phm1,phm2),
     &  'faulty internal parameter phm0, phm1, phm2')
       call ipgetr ('b2sqel_artificial_radiation', art_rad)
       call ipgeti ('b2mndr_eirene', use_eirene)
      endif
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test ev
      call xertst (0.0_R8.lt.ev, 'faulty argument ev')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
      endif
*   ..compute wrk0()=log(te()/ev)
      wrk0 = log(te/ev)

* ..main computation
*   ..compute rate coefficients
      do is = 0, ns-1
       do iCv = 1, nCv
*      ..compute ionisation rate coefficient
         rtw%rsa(iCv,is) = phm0*
     &     expu(rt%rlsa(iCv,0,is)+rt%rlsa(iCv,1,is)*wrk0(iCv))
*      ..compute recombination rate coefficient
         rtw%rra(iCv,is) = phm1*
     &     expu(rt%rlra(iCv,0,is)+rt%rlra(iCv,1,is)*wrk0(iCv))
*      ..compute heat loss rate coefficient
         rtw%rqa(iCv,is) = phm2*
     &     ev*expu(rt%rlqa(iCv,0,is)+rt%rlqa(iCv,1,is)*wrk0(iCv))
         rtw%rqa_ra(iCv,is) = phm2*
     &     ev*expu(rt%rlqa_ra(iCv,0,is)+rt%rlqa_ra(iCv,1,is)*wrk0(iCv))
*      ..compute line radiation rate coefficient
         rtw%rrd(iCv,is) = phm2*
     &     ev*expu(rt%rlrd(iCv,0,is)+rt%rlrd(iCv,1,is)*wrk0(iCv))
*      ..compute bremsstrahlung radiation rate coefficient
         rtw%rbr(iCv,is) = phm2*
     &     ev*expu(rt%rlbr(iCv,0,is)+rt%rlbr(iCv,1,is)*wrk0(iCv))
*      ..compute effective charge state
         rt%rza(iCv,is) = max(zamin(is),min(zamax(is), 
     &     rt%rlza(iCv,0,is)+rt%rlza(iCv,1,is)*wrk0(iCv)))
*      ..compute effective square charge
         rt%rz2(iCv,is) = max(zamin(is)**2,min(zamax(is)**2,
     &     rt%rlz2(iCv,0,is)+rt%rlz2(iCv,1,is)*wrk0(iCv)))
*      ..compute cumulative ionisation potential
         rt%rpt(iCv,is) = 
     &     rt%rlpt(iCv,0,is)+rt%rlpt(iCv,1,is)*wrk0(iCv)
*      ..compute effective ionisation potential
         rt%rpi(iCv,is) = 
     &     rt%rlpi(iCv,0,is)+rt%rlpi(iCv,1,is)*wrk0(iCv)
       enddo
      enddo

* ..artificial radiation                DPC 2000.08.25
      if(art_rad.ne.0.0_R8) then
        if(ncall.eq.0) then
          art_rad=abs(art_rad)
          write(*,*)
     1     'Artificial radiation model active, art_rad = ',art_rad
        endif
        is=ismain    ! XPB 2007.11.02 was ns-1
        do iCv = 1, nCv
*      ..compute artificial radiation rate coefficient
            TNORM = MAX(0.000001_R8,TE(iCv)/(15.0_R8*EV))
            rtw%rqa(iCv,is) = rtw%rqa(iCv,is) +
     1       art_rad*2.0E-31_R8/(TNORM**1.5_R8+1.0_R8/TNORM**3)
            rtw%rrd(iCv,is) = rtw%rrd(iCv,is) +
     1       art_rad*2.0E-31_R8/(TNORM**1.5_R8+1.0_R8/TNORM**3)
        enddo
      endif
      
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sqel

      end subroutine b2sqel
