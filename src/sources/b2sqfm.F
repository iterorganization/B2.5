*-----------------------------------------------------------------------
*.specification

      subroutine b2sqfm (nCv, nms, ismain, switch, ev, te, ti, rt, rtw)
      use b2mod_types
      use b2mod_math
      use b2mod_b2cmpa
      use b2us_plasma
      use b2mod_switches
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
cc      use b2mod_ad
cc     & , only : ncall_b2sqfm
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nms, ismain
      type (switches), intent (in) :: switch
      real (kind=R8) ::
     *  ev, te(nCv), ti(nCv)
      type (B2Rates), intent (inout) :: rt
*   ..output arguments (unspecified on entry)
      type (B2RatesWork), intent (inout) :: rtw

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SQFM computes on all cells of the mesh the rate coefficients
*     for fluid molecules and molecular ions.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine employs the log-log linear representation of the rate 
*     coefficients to evaluate on all cells of the mesh. 
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
      integer iCv, ims
      real (kind=R8) :: tnorm
      real (kind=R8) :: wrk0(nCv), wrk1(nCv)
*   ..procedures
      intrinsic log, min
      external xertst
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sqfm')
*   ..test nCv, nms
      call xertst (0.lt.nCv, 'faulty argument nCv')
      call xertst (1.le.nms, 'faulty argument ns')
*   ..test ev
      call xertst (0.0_R8.lt.ev, 'faulty argument ev')

*   ..compute wrk0()=log(te()/ev)
      wrk0 = log(te/ev)
*   ..compute wrk1()=log(ti()/ev) !! for elastic collisions and ion conversion
      wrk1 = log(ti/ev)

#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SQEL) )
! TNORM was originally shared in the structured code, but this seems
! wrong, as it is computed for each cell (race conditions may occur
! that way)
!$omp parallel DEFAULT(none)
!$omp& private(iCv, ims, TNORM)
!$omp& shared(wrk0, wrk1)
!$omp& shared(zamin, zamax, ev, nms, nCv)
!$omp& shared(ismain, TE)
!$omp& shared(rtw, rt, switch)
#endif

* ..main computation
*   ..compute rate coefficients
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SQEL) )
!$omp do
#endif
      do ims = 0, nms-1
       do iCv = 1, nCv
*      ..compute non-dissociative ionization rate coefficient
         rtw%r_nd(iCv,ims) = 
     &     expu(rt%rl_nd(iCv,0,ims)+rt%rl_nd(iCv,1,ims)*wrk0(iCv))

*      ..compute dissociation rate coefficient
         rtw%r_d(iCv,ims) = 
     &     expu(rt%rl_d(iCv,0,ims)+rt%rl_d(iCv,1,ims)*wrk0(iCv))

*      ..compute dissociative ionization rate coefficient
         rtw%r_di1(iCv,ims) = 
     &     expu(rt%rl_di1(iCv,0,ims)+rt%rl_di1(iCv,1,ims)*wrk0(iCv))

*      ..compute ion conversion rate coefficient
         rtw%r_ic(iCv,ims) = 
     &     expu(rt%rl_ic(iCv,0,ims)+rt%rl_ic(iCv,1,ims)*wrk1(iCv))

*      ..compute elastic collision rate coefficient (mom, parallel)
         rtw%r_ec_m(iCv,ims) = 
     &     expu(rt%rl_ec_m(iCv,0,ims)+rt%rl_ec_m(iCv,1,ims)*wrk1(iCv))

*      ..compute elastic collision rate coefficient (mom) (diamagnetic bulk velocity difference)
         rtw%r_ec_m_d(iCv,ims) = 
     &   expu(rt%rl_ec_m_d(iCv,0,ims)+rt%rl_ec_m_d(iCv,1,ims)*wrk1(iCv))

*      ..compute elastic collision rate coefficient (mom) (radial bulk velocity difference)
         rtw%r_ec_m_r(iCv,ims) = 
     &   expu(rt%rl_ec_m_r(iCv,0,ims)+rt%rl_ec_m_r(iCv,1,ims)*wrk1(iCv))

*      ..compute elastic collision rate coefficient (energy)
         rtw%r_ec_E(iCv,ims) = 
     &     expu(rt%rl_ec_E(iCv,0,ims)+rt%rl_ec_E(iCv,1,ims)*wrk1(iCv))

*      ..compute dissociative recombination rate coefficient
         rtw%r_dr(iCv,ims) = 
     &     expu(rt%rl_dr(iCv,0,ims)+rt%rl_dr(iCv,1,ims)*wrk0(iCv))

*      ..compute dissociative excitation rate coefficient
         rtw%r_de(iCv,ims) = 
     &     expu(rt%rl_de(iCv,0,ims)+rt%rl_de(iCv,1,ims)*wrk0(iCv))

*      ..compute dissociative ionization rate coefficient
         rtw%r_di2(iCv,ims) = 
     &     expu(rt%rl_di2(iCv,0,ims)+rt%rl_di2(iCv,1,ims)*wrk0(iCv))

       enddo
      enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SQEL) )
!$omp end do
#endif

#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SQEL) )
!$omp end parallel
#endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sqel

      end subroutine b2sqfm
