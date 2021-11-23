*-----------------------------------------------------------------------
*.specification

      subroutine b2spel (nCv, ns, ev, te, ne, rt)
      use b2mod_types
      use b2mod_b2cmrc
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2spel
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      real (kind=R8) ::
     *  ev, te(nCv), ne(nCv)
*   ..output arguments (unspecified on entry)
      type (B2Rates), intent (inout) :: rt
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SPEL computes on all cells of the mesh a log-log representation
*     of the rate coefficients for ionisation, recombination, radiation 
*     and electron heat loss.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine uses table interpolation to compute on each cell
*     and for each atomic species a local representation of the rate
*     coefficients for ionisation, recombination, electron heat
*     loss and radiation. The representation is linear in the log of 
*     the rate coefficient versus the log of the temperature.
*     The tables store the log of the rate coefficients as a function
*     of log(te/ev) and log(ne). Linear interpolation is used.
*     See the appropriate routines elsewhere in the code for a
*     description of the computation of the table entries.
*
*
*  8. further comments
*
*     Linear interpolation in the table is employed to represent the
*     log of each rate coefficient as a function of log(te/ev) and
*     log(ne). If log(te/ev) falls outside the domain of the table
*     then linear extrapolation in that variable is used; if log(ne)
*     falls outside the domain of the table then a prudent zeroth-order
*     extrapolation is used.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is, ite, ine
      logical ltst0
      real (kind=R8) ::
     *   rlte, rlne, fxte, fxne, t0, t1, tt
*   ..procedures
      intrinsic min, max, log
      real (kind=R8) ::
     *   trimg
      external subini, subend, xertst
      external b2xvsg, trimg
*     (trimg will be used to trim the logarithmic gradients to avoid
*     some pathetic behaviour near the underflow limit.)
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2spel')
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test ev
      call xertst (0.0_R8.lt.ev, 'faulty argument ev')
*   ..extensive tests on first few calls
      if (ncall_b2spel.lt.3) then
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
      endif

* ..compute rate coefficients
*   ..initialise indices for interpolation
      ite = 0
      ine = 0
*   ..loop over cells
! $AD II-LOOP
      do iCv = 1,nCv
*     ..set ltst0
        ltst0 = rtlt(0).le.log(te(iCv)/ev)
*     ..compute abscissae
        rlte = max(rtlt(0),min(rtlt(rtnt),log(te(iCv)/ev)))
        rlne = max(rtln(0),min(rtln(rtnn),log(ne(iCv))))
*     ..find ite
*       (linear search starting from previous ite, to establish
*         0.le.ite.lt.rtnt;
*         ite.eq.0.or.rtlt(ite).le.rlte;
*         ite.eq.rtnt-1.or.rlte.lt.rtlt(ite+1).)
    1   if (ite.lt.rtnt-1) then
         if (rtlt(ite+1).le.rlte) then
          ite = ite+1
          goto 1
         endif
        endif
    2   if (0.lt.ite) then
         if (rlte.lt.rtlt(ite)) then
          ite = ite-1
          goto 2
         endif
        endif
*     ..compute coefficient for linear interpolation
        fxte = (rlte-rtlt(ite))/(rtlt(ite+1)-rtlt(ite))
*     ..find ine
*       (linear search starting from previous ine, to establish
*         0.le.ine.lt.rtnn;
*         ine.eq.0.or.rtln(ine).le.rlne;
*         ine.eq.rtnn-1.or.rlne.lt.rtln(ine+1).)
    3   if (ine.lt.rtnn-1) then
         if (rtln(ine+1).le.rlne) then
          ine = ine+1
          goto 3
         endif
        endif
    4   if (0.lt.ine) then
         if (rlne.lt.rtln(ine)) then
          ine = ine-1
          goto 4
         endif
        endif
*     ..compute coefficient for linear interpolation
        fxne = (rlne-rtln(ine))/(rtln(ine+1)-rtln(ine))
*     ..compute log-log rate coefficients
        do is = 0, ns-1
*      ..compute ionisation rate coefficient
         t0 = (1.0_R8-fxne)*rtlsa(ite,ine,is)+
     &        fxne*rtlsa(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlsa(ite+1,ine,is)+
     &        fxne*rtlsa(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rt%rlsa(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rt%rlsa(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rt%rlsa(iCv,0,is) = 
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rlsa(iCv,1,is)*rlte
*      ..compute recombination rate coefficient
         t0 = (1.0_R8-fxne)*rtlra(ite,ine,is)+
     &        fxne*rtlra(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlra(ite+1,ine,is)+
     &        fxne*rtlra(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rt%rlra(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rt%rlra(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rt%rlra(iCv,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rlra(iCv,1,is)*rlte
*      ..compute heat loss rate coefficient
         t0 = (1.0_R8-fxne)*rtlqa(ite,ine,is)+
     &        fxne*rtlqa(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlqa(ite+1,ine,is)+
     &        fxne*rtlqa(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rt%rlqa(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rt%rlqa(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rt%rlqa(iCv,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rlqa(iCv,1,is)*rlte 
         t0 = (1.0_R8-fxne)*rtlqa_ra(ite,ine,is)+
     &        fxne*rtlqa_ra(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlqa_ra(ite+1,ine,is)+
     &        fxne*rtlqa_ra(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rt%rlqa_ra(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rt%rlqa_ra(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rt%rlqa_ra(iCv,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rlqa_ra(iCv,1,is)*rlte
*      ..compute line radiation rate coefficient
         t0 = (1.0_R8-fxne)*rtlrd(ite,ine,is)+
     &        fxne*rtlrd(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlrd(ite+1,ine,is)+
     &        fxne*rtlrd(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rt%rlrd(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rt%rlrd(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rt%rlrd(iCv,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rlrd(iCv,1,is)*rlte
*      ..compute bremsstrahlung radiation rate coefficient
         t0 = (1.0_R8-fxne)*rtlbr(ite,ine,is)+
     &        fxne*rtlbr(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlbr(ite+1,ine,is)+
     &        fxne*rtlbr(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rt%rlbr(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rt%rlbr(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rt%rlbr(iCv,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rlbr(iCv,1,is)*rlte
*      ..compute effective charge state
         t0 = (1.0_R8-fxne)*rtlza(ite,ine,is)+
     &        fxne*rtlza(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlza(ite+1,ine,is)+
     &        fxne*rtlza(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rt%rlza(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rt%rlza(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rt%rlza(iCv,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rlza(iCv,1,is)*rlte
*      ..compute effective square charge
         t0 = (1.0_R8-fxne)*rtlz2(ite,ine,is)+
     &        fxne*rtlz2(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlz2(ite+1,ine,is)+
     &        fxne*rtlz2(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rt%rlz2(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rt%rlz2(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rt%rlz2(iCv,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rlz2(iCv,1,is)*rlte
*      ..compute cumulative ionisation potential
         t0 = (1.0_R8-fxne)*rtlpt(ite,ine,is)+
     &        fxne*rtlpt(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlpt(ite+1,ine,is)+
     &        fxne*rtlpt(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rt%rlpt(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rt%rlpt(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rt%rlpt(iCv,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rlpt(iCv,1,is)*rlte
*      ..compute effective ionisation potential
         t0 = (1.0_R8-fxne)*rtlpi(ite,ine,is)+
     &        fxne*rtlpi(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlpi(ite+1,ine,is)+
     &        fxne*rtlpi(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rt%rlpi(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rt%rlpi(iCv,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rt%rlpi(iCv,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rlpi(iCv,1,is)*rlte
        enddo
      enddo

* ..return
      ncall_b2spel = ncall_b2spel+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Note: This routine should be rewritten for vectorisation.

*-----------------------------------------------------------------------
*.end b2spel

      end subroutine b2spel

