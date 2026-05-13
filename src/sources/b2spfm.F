*-----------------------------------------------------------------------
*.specification

      subroutine b2spfm (nCv, nms, te, ti, ne, ua, rt, dv)
      use b2mod_types
      use b2mod_b2cmrc_mol
      use b2us_plasma
      use b2mod_math
      use b2mod_openmp
      use b2mod_subsys
      use b2mod_fluid_molecules
      implicit none
*   ..input arguments (unchanged on exit)
      type (B2Derivatives), intent (inout) :: dv
      integer nCv, nms
      real (kind=R8) ::
     *  te(nCv), ne(nCv), ti(nCv), ua(nCv)
*   ..output arguments (unspecified on entry)
      type (B2Rates), intent (inout) :: rt
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SPFM computes on all cells of the mesh a log-log representation
*     of the rate coefficients for neutral molecules and molecular ions.
*
*
*  3. description
*
*     This routine uses table interpolation to compute on each cell
*     and for each molecular species a local representation of the rate
*     coefficients. The representation is linear in the log of
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
      integer iCv, ims, ite, ine, itr, iti, 
     *        i_mcoll_p, i_mcoll_d, i_mcoll_r
      logical ltst0
      real (kind=R8) ::
     *   rlte, rlne, fxte, fxne, fxtr, t0, t1, tt, vth, tr,
     *   rltr, rlti, fxti, wrk(0:rtn_mcoll),
     *   mcoll_p, mcoll_d, mcoll_r, rl_mcoll_p, rl_mcoll_d, rl_mcoll_r
*   ..procedures
      intrinsic min, max, log
      external xertst
      external b2xvsg
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2spfm')
*   ..test nCv, ns
      call xertst (0.lt.nCv, 'faulty argument nCv')
      call xertst (1.le.nms, 'faulty argument nms')
*   ..test ev
      call xertst (0.0_R8.lt.ev, 'faulty argument ev')

      if (background_was_read.eq.0) then
        call read_molecules_background_ext(nCv)
      endif

* ..compute rate coefficients
*   ..initialise indices for interpolation
      ite = 0
      ine = 0
      iti = 0
*   ..loop over cells
! $AD II-LOOP
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SPEL) )
!$omp parallel do default(none) shared(nCv,nms) private(iCv,ims)
!$omp& firstprivate(ite, ine, iti)
!$omp& private(rlte, rlne, fxte, fxne, t0, t1, tt, fxti, rlti)
!$omp& shared(rtlt, te, ne, ti, rtnt, rtln, rtnn)
!$omp& shared(rtl_nd, rtl_d, rtl_di1, rtl_ic, rtl_ec_m, rtl_ec_E)
!$omp& shared(rtl_dr, rtl_de, rtl_di2)
!$omp& shared(rt)
#endif
      do iCv = 1,nCv
*     ..compute abscissae
        rlte = max(rtlt(0),min(rtlt(rtnt),log(te(iCv)/ev)))
        rlne = max(rtln(0),min(rtln(rtnn),log(ne(iCv))))
        rlti = max(rtlt(0),min(rtlt(rtnt),log(ti(iCv)/ev)))
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

*     ..find iti (only needed for the ion conversion reaction, which depends on Ti instead of Te)
*       (linear search starting from previous it, to establish
*         0.le.iti.lt.rtnt;
*         iti.eq.0.or.rtlt(iti).le.rlti;
*         iti.eq.rtnt-1.or.rlti.lt.rtlt(iti+1).)
   11   if (iti.lt.rtnt-1) then
         if (rtlt(iti+1).le.rlti) then
          iti = iti+1
          goto 11
         endif
        endif
   22   if (0.lt.iti) then
         if (rlti.lt.rtlt(iti)) then
          iti = iti-1
          goto 22
         endif
        endif
*     ..compute coefficient for linear interpolation
        fxti = (rlti-rtlt(iti))/(rtlt(iti+1)-rtlt(iti))

*     ..compute log-log rate coefficients
        do ims = 0, nms-1
*      ..compute non-dissociative ionization rate coefficient
         t0 = (1.0_R8-fxne)*rtl_nd(ite,ine,ims)+
     &        fxne*rtl_nd(ite,ine+1,ims)
         t1 = (1.0_R8-fxne)*rtl_nd(ite+1,ine,ims)+
     &        fxne*rtl_nd(ite+1,ine+1,ims)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         rt%rl_nd(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt)
         rt%rl_nd(iCv,0,ims) =
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rl_nd(iCv,1,ims)*rlte

*      ..compute dissociation rate coefficient
         t0 = (1.0_R8-fxne)*rtl_d(ite,ine,ims)+
     &        fxne*rtl_d(ite,ine+1,ims)
         t1 = (1.0_R8-fxne)*rtl_d(ite+1,ine,ims)+
     &        fxne*rtl_d(ite+1,ine+1,ims)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         rt%rl_d(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt)
         rt%rl_d(iCv,0,ims) =
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rl_d(iCv,1,ims)*rlte

*      ..compute dissociative ionization rate coefficient
         t0 = (1.0_R8-fxne)*rtl_di1(ite,ine,ims)+
     &        fxne*rtl_di1(ite,ine+1,ims)
         t1 = (1.0_R8-fxne)*rtl_di1(ite+1,ine,ims)+
     &        fxne*rtl_di1(ite+1,ine+1,ims)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         rt%rl_di1(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt)
         rt%rl_di1(iCv,0,ims) =
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rl_di1(iCv,1,ims)*rlte

*      ..compute ion conversion rate coefficient !! iti instead of ite
         t0 = (1.0_R8-fxne)*rtl_ic(iti,ine,ims)+
     &        fxne*rtl_ic(iti,ine+1,ims)
         t1 = (1.0_R8-fxne)*rtl_ic(iti+1,ine,ims)+
     &        fxne*rtl_ic(iti+1,ine+1,ims)
         tt = (t1-t0)/(rtlt(iti+1)-rtlt(iti))
         rt%rl_ic(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt) 
         rt%rl_ic(iCv,0,ims) =
     .    (1.0_R8-fxti)*t0+fxti*t1-rt%rl_ic(iCv,1,ims)*rlti

*      ..compute dissociative recombination
         t0 = (1.0_R8-fxne)*rtl_dr(ite,ine,ims)+
     &        fxne*rtl_dr(ite,ine+1,ims)
         t1 = (1.0_R8-fxne)*rtl_dr(ite+1,ine,ims)+
     &        fxne*rtl_dr(ite+1,ine+1,ims)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         rt%rl_dr(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt)
         rt%rl_dr(iCv,0,ims) =
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rl_dr(iCv,1,ims)*rlte

*      ..compute dissociative excitation
         t0 = (1.0_R8-fxne)*rtl_de(ite,ine,ims)+
     &        fxne*rtl_de(ite,ine+1,ims)
         t1 = (1.0_R8-fxne)*rtl_de(ite+1,ine,ims)+
     &        fxne*rtl_de(ite+1,ine+1,ims)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         rt%rl_de(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt)
         rt%rl_de(iCv,0,ims) =
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rl_de(iCv,1,ims)*rlte

*      ..compute dissociative ionization (from molecular ion)
         t0 = (1.0_R8-fxne)*rtl_di2(ite,ine,ims)+
     &        fxne*rtl_di2(ite,ine+1,ims)
         t1 = (1.0_R8-fxne)*rtl_di2(ite+1,ine,ims)+
     &        fxne*rtl_di2(ite+1,ine+1,ims)
         tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
         rt%rl_di2(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt)
         rt%rl_di2(iCv,0,ims) =
     .    (1.0_R8-fxte)*t0+fxte*t1-rt%rl_di2(iCv,1,ims)*rlte

        enddo !! end loop over molecular species
      enddo !! end loop over cells
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SPEL) )
!$omp end parallel do
#endif

cc more complex treatment for elastic collisions
*   ..initialise indices for interpolation
      ite = 0
      itr = 0
*   ..loop over cells
      do iCv = 1,nCv
*     ..compute abscissae
        tr = t_mol(iCv,0)/ti(iCv) !! ratio of molecule to ion temperature
        vth = sqrt(2.0_R8*ti(iCv)/(mp*am(1))) !! TODO generalize species
        mcoll_p = abs(u_mol(iCv,0)-ua(iCv))/vth !! velocity difference in parallel direction expressed as Mach number   
        mcoll_d = abs(u_pol_mol(iCv,0))/vth !! velocity difference in diamagnetic direction expressed as Mach number  (neglecting ion component)
        mcoll_r = abs(u_rad_mol(iCv,0))/vth !! velocity difference in radial direction expressed as Mach number  (neglecting ion component)
        rlti = max(rtlt(0),min(rtlt(rtnt),log(ti(iCv)/ev)))
        rltr = max(rtltr(0),min(rtltr(rtntr),log(tr)))
        rl_mcoll_p = max(rtl_mcoll(0),min(rtl_mcoll(rtn_mcoll),
     .              log(mcoll_p)))
        rl_mcoll_d = max(rtl_mcoll(0),min(rtl_mcoll(rtn_mcoll),
     .              log(mcoll_d)))
        rl_mcoll_r = max(rtl_mcoll(0),min(rtl_mcoll(rtn_mcoll),
     .              log(mcoll_r)))
*     ..find iti
*       (linear search starting from previous iti, to establish
*         0.le.ite.lt.rtnt;
*         iti.eq.0.or.rtlt(iti).le.rlti;
*         iti.eq.rtnt-1.or.rlti.lt.rtlt(iti+1).)
    5   if (iti.lt.rtnt-1) then
         if (rtlt(iti+1).le.rlti) then
          iti = iti+1
          goto 5
         endif
        endif
    6   if (0.lt.iti) then
         if (rlti.lt.rtlt(iti)) then
          iti = iti-1
          goto 6
         endif
        endif
*     ..compute coefficient for linear interpolation
        fxti = (rlti-rtlt(iti))/(rtlt(iti+1)-rtlt(iti))
*     ..find itr
*       (linear search starting from previous itr, to establish
*         0.le.itr.lt.rtntr;
*         itr.eq.0.or.rtltr(itr).le.rltr;
*         itr.eq.rtntr-1.or.rltr.lt.rtltr(itr+1).)
    7   if (itr.lt.rtntr-1) then
         if (rtltr(itr+1).le.rltr) then
          itr = itr+1
          goto 7
         endif
        endif
    8   if (0.lt.itr) then
         if (rltr.lt.rtltr(itr)) then
          itr = itr-1
          goto 8
         endif
        endif
*     ..compute coefficient for linear interpolation
        fxtr = (rltr-rtltr(itr))/(rtltr(itr+1)-rtltr(itr))
cc for the velocity difference, we use nearest neighbour instead of interpolation
*     ..find i_mcoll
cc parallel
       wrk = rtl_mcoll-rl_mcoll_p
       wrk = abs(wrk)
       i_mcoll_p = MINLOC(wrk,1)-1 !! MINLOC seems to ignore the 0-indexing (returning 1 if the first element is the smallest, even though that has index 0) -> -1
cc diamagnetic
       wrk = rtl_mcoll-rl_mcoll_d
       wrk = abs(wrk)
       i_mcoll_d = MINLOC(wrk,1)-1 !! MINLOC seems to ignore the 0-indexing (returning 1 if the first element is the smallest, even though that has index 0) -> -1
cc radial
       wrk = rtl_mcoll-rl_mcoll_r
       wrk = abs(wrk)
       i_mcoll_r = MINLOC(wrk,1)-1 !! MINLOC seems to ignore the 0-indexing (returning 1 if the first element is the smallest, even though that has index 0) -> -1

*     ..compute log-log rate coefficients
        do ims = 0, nms-1

*      ..compute elastic collision rate coefficient (parallel momentum)
         t0 = (1.0_R8-fxtr)*rtl_ec_m(iti,itr,i_mcoll_p,ims)+
     &        fxtr*rtl_ec_m(iti,itr+1,i_mcoll_p,ims)
         t1 = (1.0_R8-fxtr)*rtl_ec_m(iti+1,itr,i_mcoll_p,ims)+
     &        fxtr*rtl_ec_m(iti+1,itr+1,i_mcoll_p,ims)
         tt = (t1-t0)/(rtlt(iti+1)-rtlt(iti))
         rt%rl_ec_m(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt)
         rt%rl_ec_m(iCv,0,ims) =
     &    (1.0_R8-fxti)*t0+fxti*t1-rt%rl_ec_m(iCv,1,ims)*rlti

*      ..compute elastic collision rate coefficient (diamagnetic momentum)
         t0 = (1.0_R8-fxtr)*rtl_ec_m(iti,itr,i_mcoll_d,ims)+
     &        fxtr*rtl_ec_m(iti,itr+1,i_mcoll_d,ims)
         t1 = (1.0_R8-fxtr)*rtl_ec_m(iti+1,itr,i_mcoll_d,ims)+
     &        fxtr*rtl_ec_m(iti+1,itr+1,i_mcoll_d,ims)
         tt = (t1-t0)/(rtlt(iti+1)-rtlt(iti))
         rt%rl_ec_m_d(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt)
         rt%rl_ec_m_d(iCv,0,ims) =
     &    (1.0_R8-fxti)*t0+fxti*t1-rt%rl_ec_m_d(iCv,1,ims)*rlti

*      ..compute elastic collision rate coefficient (radial momentum)
         t0 = (1.0_R8-fxtr)*rtl_ec_m(iti,itr,i_mcoll_r,ims)+
     &        fxtr*rtl_ec_m(iti,itr+1,i_mcoll_r,ims)
         t1 = (1.0_R8-fxtr)*rtl_ec_m(iti+1,itr,i_mcoll_r,ims)+
     &        fxtr*rtl_ec_m(iti+1,itr+1,i_mcoll_r,ims)
         tt = (t1-t0)/(rtlt(iti+1)-rtlt(iti))
         rt%rl_ec_m_r(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt)
         rt%rl_ec_m_r(iCv,0,ims) =
     &    (1.0_R8-fxti)*t0+fxti*t1-rt%rl_ec_m_r(iCv,1,ims)*rlti

*      ..compute elastic collision rate coefficient (energy)
         t0 = (1.0_R8-fxtr)*rtl_ec_E(iti,itr,ims)+
     &        fxtr*rtl_ec_E(iti,itr+1,ims)
         t1 = (1.0_R8-fxtr)*rtl_ec_E(iti+1,itr,ims)+
     &        fxtr*rtl_ec_E(iti+1,itr+1,ims)
         tt = (t1-t0)/(rtlt(iti+1)-rtlt(iti))
         rt%rl_ec_E(iCv,1,ims) = trimg(-20.0_R8,20.0_R8,tt)
         rt%rl_ec_E(iCv,0,ims) =
     .    (1.0_R8-fxti)*t0+fxti*t1-rt%rl_ec_E(iCv,1,ims)*rlti


        enddo
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Note: This routine should be rewritten for vectorisation.

*-----------------------------------------------------------------------
*.end b2spel

      end subroutine b2spfm
