*-----------------------------------------------------------------------
*.specification

      subroutine b2spel (nx, ny, ns, ev, te, ne, 
     &  rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi)
      use b2mod_types
      use b2mod_b2cmrc
#ifdef WITH_FTIMINGS
      use timings
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  ev, te(-1:nx,-1:ny), ne(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  rlsa(-1:nx,-1:ny,0:1,0:ns-1), rlra(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlqa(-1:nx,-1:ny,0:1,0:ns-1), rlrd(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlbr(-1:nx,-1:ny,0:1,0:ns-1), rlza(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlz2(-1:nx,-1:ny,0:1,0:ns-1), rlpt(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlpi(-1:nx,-1:ny,0:1,0:ns-1)
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
      integer ncall, n2, ix, iy, is, ite, ine
      logical ltst0
      save ncall
      real (kind=R8) ::
     *   rlte, rlne, fxte, fxne, t0, t1, tt
*   ..procedures
      intrinsic min, max, log
      real (kind=R8) ::
     *   trimg
      external xertst
      external b2xvsg
      trimg(t0,t1,tt) = max(t0,min(t1,tt))
*     (trimg will be used to trim the logarithmic gradients to avoid
*     some pathetic behaviour near the underflow limit.)
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2spel')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test ev
      call xertst (0.0_R8.lt.ev, 'faulty argument ev')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      endif

* ..compute rate coefficients
*   ..initialise indices for interpolation
      ite = 0
      ine = 0
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
*   ..loop over cells
#ifndef NO_OPENMP_B2SPEL
!$OMP PARALLEL DO default(none) shared(nx,ny,ns) private(ix,iy,is)
!$OMP& FIRSTPRIVATE(ite, ine)
!$OMP& private(ltst0, rlte, rlne, fxte, fxne, t0, t1, tt)
!$OMP& shared(rtlt, te, ev, rtnt, rtln, rtnn, rtlsa, rlsa, rtlra, rlra)
!$OMP& shared(ne, rtlqa, rlqa, rtlrd, rlrd, rtlbr, rlbr, rtlza, rlza)
!$OMP& shared(rtlz2, rlz2, rtlpt, rlpt, rtlpi, rlpi)
#endif
      do is = 0, ns-1
*     The IS loop is moved outside, so that we can calculate parallel in
*     IS. This will ensure that the output arrays are distributed along
*     the slowest varying index is.
      do iy = -1, ny
       do ix = -1, nx
*     ..set ltst0
        ltst0 = rtlt(0).le.log(te(ix,iy)/ev)
*     ..compute abscissae
        rlte = max(rtlt(0),min(rtlt(rtnt),log(te(ix,iy)/ev)))
        rlne = max(rtln(0),min(rtln(rtnn),log(ne(ix,iy))))
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
*      ..compute ionisation rate coefficient
         t0 = (1.0_R8-fxne)*rtlsa(ite,ine,is)+
     &        fxne*rtlsa(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlsa(ite+1,ine,is)+
     &        fxne*rtlsa(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rlsa(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rlsa(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rlsa(ix,iy,0,is) = 
     .    (1.0_R8-fxte)*t0+fxte*t1-rlsa(ix,iy,1,is)*rlte
*      ..compute recombination rate coefficient
         t0 = (1.0_R8-fxne)*rtlra(ite,ine,is)+
     &        fxne*rtlra(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlra(ite+1,ine,is)+
     &        fxne*rtlra(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rlra(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rlra(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rlra(ix,iy,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rlra(ix,iy,1,is)*rlte
*      ..compute heat loss rate coefficient
         t0 = (1.0_R8-fxne)*rtlqa(ite,ine,is)+
     &        fxne*rtlqa(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlqa(ite+1,ine,is)+
     &        fxne*rtlqa(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rlqa(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rlqa(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rlqa(ix,iy,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rlqa(ix,iy,1,is)*rlte
*      ..compute line radiation rate coefficient
         t0 = (1.0_R8-fxne)*rtlrd(ite,ine,is)+
     &        fxne*rtlrd(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlrd(ite+1,ine,is)+
     &        fxne*rtlrd(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rlrd(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rlrd(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rlrd(ix,iy,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rlrd(ix,iy,1,is)*rlte
*      ..compute bremsstrahlung radiation rate coefficient
         t0 = (1.0_R8-fxne)*rtlbr(ite,ine,is)+
     &        fxne*rtlbr(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlbr(ite+1,ine,is)+
     &        fxne*rtlbr(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rlbr(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rlbr(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rlbr(ix,iy,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rlbr(ix,iy,1,is)*rlte
*      ..compute effective charge state
         t0 = (1.0_R8-fxne)*rtlza(ite,ine,is)+
     &        fxne*rtlza(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlza(ite+1,ine,is)+
     &        fxne*rtlza(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rlza(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rlza(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rlza(ix,iy,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rlza(ix,iy,1,is)*rlte
*      ..compute effective square charge
         t0 = (1.0_R8-fxne)*rtlz2(ite,ine,is)+
     &        fxne*rtlz2(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlz2(ite+1,ine,is)+
     &        fxne*rtlz2(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rlz2(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rlz2(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rlz2(ix,iy,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rlz2(ix,iy,1,is)*rlte
*      ..compute cumulative ionisation potential
         t0 = (1.0_R8-fxne)*rtlpt(ite,ine,is)+
     &        fxne*rtlpt(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlpt(ite+1,ine,is)+
     &        fxne*rtlpt(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rlpt(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rlpt(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rlpt(ix,iy,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rlpt(ix,iy,1,is)*rlte
*      ..compute effective ionisation potential
         t0 = (1.0_R8-fxne)*rtlpi(ite,ine,is)+
     &        fxne*rtlpi(ite,ine+1,is)
         t1 = (1.0_R8-fxne)*rtlpi(ite+1,ine,is)+
     &        fxne*rtlpi(ite+1,ine+1,is)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         if (ltst0) then
          rlpi(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         else
          rlpi(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         endif
         rlpi(ix,iy,0,is) =  
     .    (1.0_R8-fxte)*t0+fxte*t1-rlpi(ix,iy,1,is)*rlte
        enddo ! ix = -1, nx
       enddo ! iy = -1, ny
      enddo ! is = 0, ns - 1
#ifndef NO_OPENMP_B2SPEL
!$OMP END PARALLEL DO
#endif
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Note: This routine should be rewritten for vectorisation.

*-----------------------------------------------------------------------
*.end b2spel

      end subroutine b2spel

