*-----------------------------------------------------------------------
*.specification

      subroutine b2spcx (nx, ny, ns, ev, amh0, th0, ne, rlcx0)
      use b2mod_types
      use b2mod_math
      use b2mod_b2cmrc
#ifdef WITH_FTIMINGS
      use timings
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *   ev, amh0, th0(-1:nx,-1:ny), ne(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   rlcx0(-1:nx,-1:ny,0:1,0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SPCX computes on all cells of the mesh a log-log representation
*     of the rate coefficients for charge exchange.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine uses table interpolation to compute on each cell
*     and for each atomic species a local representation of the rate
*     coefficient for charge exchange with neutral hydrogen.
*     The representation is linear in the log of the rate coefficient
*     versus the log of the atom temperature per amu.
*     The table stores the log of the rate coefficient as a function of
*     the log of the neutral hydrogen temperature per atomic mass and
*     the log of the electron density. Linear interpolation is used.
*     See the appropriate routines elsewhere in the code for a
*     description of the computation of the table entries.
*
*
*  8. further comments
*
*     Linear interpolation in the table is employed to represent the
*     log of the rate coefficient as a function of log(th0/(amh0*ev0))
*     and log(nh0). If log(th0/(amh0*ev)) falls outside the domain of
*     the table then linear extrapolation in that variable is used;
*     if log(nh0) falls outside the domain of the table then a prudent
*     zeroth-order extrapolation is used.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, it0, in0
      save ncall
      real (kind=R8) ::
     *   rlt0, rln0, fxt0, fxn0, t0, t1, tt
*   ..procedures
      intrinsic log
      external subini, subend, xertst
      external b2xvsg
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2spcx')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test ev, amh0
      call xertst (0.lt.ev, 'faulty argument ev')
      call xertst (0.lt.amh0.and.amh0.le.3,
     &  'faulty argument range amh0')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
       call b2xvsg (n2, th0, 1, 'th0', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      endif

* ..compute rate coefficients
*   ..initialise indices for interpolation
      it0 = 0
      in0 = 0
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
*   ..loop over cells
#ifndef NO_OPENMP_B2SPCX
!$OMP PARALLEL DO
!$OMP& default(none)
!$OMP& firstprivate(it0, in0)
!$OMP& private(ix, iy, is)
!$OMP& private(rlt0, rln0, fxt0, fxn0, t0, t1, tt)
!$OMP& shared(rtlt, th0)
!$OMP& shared(nx, ny, ns)
!$OMP& shared(amh0, ev, rtnt, rtln, rtnn, ne, rtlcx)
! this is the result variable of the loop:
!$OMP& shared(rlcx0)
#endif
      do is = 0, ns-1
       do iy = -1, ny
        do ix = -1, nx
*      ..compute abscissae
         rlt0 = max(rtlt(0),min(rtlt(rtnt),log(th0(ix,iy)/(amh0*ev))))
         rln0 = max(rtln(0),min(rtln(rtnn),log(ne(ix,iy))))
*      ..find it0
*        (linear search starting from previous it0, to establish
*          0.le.it0.lt.rtnt;
*          it0.eq.0.or.rtlt(it0).le.rlt0;
*          it0.eq.rtnt-1.or.rlt0.lt.rtlt(it0+1).)
    1    if (it0.lt.rtnt-1) then
          if (rtlt(it0+1).le.rlt0) then
           it0 = it0+1
           goto 1
          endif
         endif
    2    if (0.lt.it0) then
          if (rlt0.lt.rtlt(it0)) then
           it0 = it0-1
           goto 2
          endif
         endif
*      ..compute coefficient for linear interpolation
         fxt0 = (rlt0-rtlt(it0))/(rtlt(it0+1)-rtlt(it0))
*      ..find in0
*        (linear search starting from previous in0, to establish
*          0.le.in0.lt.rtnn;
*          in0.eq.0.or.rtln(in0).le.rln0;
*          in0.eq.rtnn-1.or.rln0.lt.rtln(in0+1).)
    3    if (in0.lt.rtnn-1) then
          if (rtln(in0+1).le.rln0) then
           in0 = in0+1
           goto 3
          endif
         endif
    4    if (0.lt.in0) then
          if (rln0.lt.rtln(in0)) then
           in0 = in0-1
           goto 4
          endif
         endif
*      ..compute coefficient for linear interpolation
         fxn0 = (rln0-rtln(in0))/(rtln(in0+1)-rtln(in0))
*      ..compute log-log rate coefficients
*      ..compute charge exchange rate coefficient
         t0 = (1-fxn0)*rtlcx(it0,in0,is)+fxn0*rtlcx(it0,in0+1,is)
         t1 = (1-fxn0)*rtlcx(it0+1,in0,is)+fxn0*rtlcx(it0+1,in0+1,is)
         tt = (t1-t0)/(rtlt(it0+1)-rtlt(it0))
         rlcx0(ix,iy,1,is) = trimg(-20.0_R8,20.0_R8,tt)
         rlcx0(ix,iy,0,is) = (1-fxt0)*t0+fxt0*t1-rlcx0(ix,iy,1,is)*rlt0
        enddo
       enddo
      enddo
#ifndef NO_OPENMP_B2SPCX
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
*.end b2spcx

      end subroutine b2spcx
