*-----------------------------------------------------------------------
*.specification

      subroutine b2upco (nCv, nFc, nVx, geo, mpg, switch, nregionv, 
     &   ua_solve, solvereg, rxf, pcm0, corpb, pccb0, pb, nb, ub, isb)
      use b2mod_types
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2mod_switches
      use b2mod_numerics_namelist
     , , only : corr_core_dn
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2upco
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, nregionv, isb
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
      real (kind=R8) ::
     *  rxf, pcm0, corpb(nCv), pccb0(nCv)
      logical ua_solve, solvereg(0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *  pb(nCv), nb(nCv), ub(nCv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UPCO corrects the state of the plasma towards the solution of
*     continuity equation for one atomic species.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer i, iCv, iFc, iFt, ireg
*     ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(nFc), wrk(nCv), wrkv(nVx)
      real (kind=R8) :: t0,
     *  pb_(nCv), nb_(nCv) ,                                        !lk 30.06.17
     *  tt0(mpg%nFt), tt1(mpg%nFt), tt2(mpg%nFt)                    !lk 30.06.17
      logical too_low_density                                       ! IYS 27.11.2017
*   ..procedures
      external subini, subend, xertst, ipgetr, b2xxgs
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2upco')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..test rxf, pcm0
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     .  'faulty argument rxf')
      call xertst (0.0_R8.le.pcm0, 'faulty argument pcm0')

* ..implement the correction
*   ..compute wrk
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, corpb, wrkv, wrk)
*   ..compute new pb and nb
      do iCv = 1, nCv
        if (0.0_R8.le.corpb(iCv)) then
         t0 = 1.0_R8+rxf*corpb(iCv)/pb(iCv)
        else
         t0 = 1.0_R8/(1.0_R8-rxf*corpb(iCv)/pb(iCv))
        endif
        nb_(iCv) = t0*nb(iCv)                        !lk 30.06.2017
        pb_(iCv) = t0*pb(iCv)                        !lk 30.06.2017
      enddo
      if(corr_core_dn(isb).ne.1.0_R8 .and.
     &               (.not. is_neutral(isb))) then  
*   ..smooth the correction for inner flux surfaces  !lk 30.06.2017{
        tt0 = 0.0e0_R8
        tt1 = 0.0e0_R8
        tt2 = 0.0e0_R8
        do iFt = 1, mpg%nFt
          do i = mpg%ftCvP(iFt,1), mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
            iCv = mpg%ftCv(i)
            if (mpg%cvOnClosedSurface(iCv).and.(iCv.le.mpg%nCi)) then ! IYS 27.03.2019
              tt0(iFt) = tt0(iFt) + geo%cvVol(iCv)
              tt1(iFt) = tt1(iFt) + geo%cvVol(iCv)*(nb_(iCv)-nb(iCv))
              tt2(iFt) = tt2(iFt) + geo%cvVol(iCv)*(pb_(iCv)-pb(iCv))
            endif
           enddo
        enddo
        do iFt = 1, mpg%nFt
          too_low_density = .FALSE.     ! IYS 27.11.2017
! Check whether the density and/or pressure falls below the threshold 'b2mndr_na_min'
! If it happens, cancel the smoothing of the correction for the whole
! flux surface where such event happens.
          do i = mpg%ftCvP(iFt,1), mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
            iCv = mpg%ftCv(i)
            if (mpg%cvOnClosedSurface(iCv).and.(iCv.le.mpg%nCi)) then ! IYS 27.03.2019
              t0 = pb(iCv)+tt2(iFt)/tt0(iFt) +
     *         corr_core_dn(isb)*(pb_(iCv)-pb(iCv)-tt2(iFt)/tt0(iFt)) ! IYS 19.12.2017
              if (t0 .lt. pb(iCv)*switch%b2mndr_na_min / nb(iCv)) then
                too_low_density = .TRUE.
                exit
              else
                pb(iCv) = t0
              endif
              t0 = nb(iCv)+tt1(iFt)/tt0(iFt) +
     *         corr_core_dn(isb)*(nb_(iCv)-nb(iCv)-tt1(iFt)/tt0(iFt))   ! IYS 19.12.2017
              if (t0 .lt. switch%b2mndr_na_min) then
                too_low_density = .TRUE.
                exit
              else
                nb(iCv) = t0
              endif
            else                                      !srv 25.07.17 {
              nb(iCv) = nb_(iCv)
              pb(iCv) = pb_(iCv)                      !srv 25.07.17 }
            endif
          enddo
          if (too_low_density) then  ! cancel smoothing of the correction on the surface if the corrected density is too low
            do i = mpg%ftCvP(iFt,1), mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
              iCv = mpg%ftCv(i)
              nb(iCv) = nb_(iCv)
              pb(iCv) = pb_(iCv)                      !srv 25.07.17 }
            enddo
          endif
        enddo
      else         ! do not smooth the correction
        nb = nb_
        pb = pb_
      endif                                          !srv 25.07.17 }                                           !srv 25.07.17
*   ..compute new ub
      if (ua_solve) then
       do iCv = 1, mpg%nCi ! internal cells only
*       (note the further damping on the velocity correction)
         ireg = mpg%cvReg(iCv)
         if(solvereg(ireg))
     &    ub(iCv) = ub(iCv)-pcm0*rxf*pccb0(iCv)*wrk(iCv)
       enddo
      endif

* ..return
      ncall_b2upco = ncall_b2upco+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2upco

      end
