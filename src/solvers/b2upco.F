*-----------------------------------------------------------------------
*.specification

      subroutine b2upco (nCv, nFc, nVx, geo, mpg, switch, nregionv, 
     &   ua_solve, solvereg, rxf, pcm0, corpb, pccb0, pb, nb, ub)
      use b2mod_types
      use b2us_geo
      use b2us_map
      use b2mod_switches
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2upco
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, nregionv
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
      integer iCv, iFc, ireg
*     ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(nFc), wrk(nCv), wrkv(nVx)
      real (kind=R8) :: t0,
     *  pb_(nCv), nb_(nCv) !,                                               !lk 30.06.17
!WG_TODO     *  tt0(-1:ny), tt1(-1:ny), tt2(-1:ny)                                !lk 30.06.17
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
      if(switch%c_corr_core_dn.eq.1.0_R8) then                  !srv 25.07.17 {
        nb = nb_
        pb = pb_
      else                                               !srv 25.07.17 }
        call xerrab('c_corr_core_dn .ne. 1 not yet available for WG')
#ifdef WG_TODO
*   ..smooth the correction for inner flux surfaces      !lk 30.06.2017{
        tt0 = 0.0e0_R8
        tt1 = 0.0e0_R8
        tt2 = 0.0e0_R8
        do iy = -1,ny
          do ix = -1,nx
            if(on_closed_surface(ix,iy).and.bottomiy(ix,iy).ne.-2) then
              tt0(iy) = tt0(iy) + vol(ix,iy)
              tt1(iy) = tt1(iy) + vol(ix,iy)*(nb_(ix,iy)-nb(ix,iy))
              tt2(iy) = tt2(iy) + vol(ix,iy)*(pb_(ix,iy)-pb(ix,iy))
            endif
          enddo
        enddo
        do iy = -1,ny
          do ix = -1,nx
            if(on_closed_surface(ix,iy).and.bottomiy(ix,iy).ne.-2) then
              nb(ix,iy) = nb(ix,iy)+tt1(iy)/tt0(iy) +
     *         switch%c_corr_core_dn*(nb_(ix,iy)-nb(ix,iy)-tt1(iy)/tt0(iy))
              pb(ix,iy) = pb(ix,iy)+tt2(iy)/tt0(iy) +
     *         switch%c_corr_core_dn*(pb_(ix,iy)-pb(ix,iy)-tt2(iy)/tt0(iy))
            else                                          !srv 25.07.17 {
              nb(ix,iy) = nb_(ix,iy)
              pb(ix,iy) = pb_(ix,iy)                      !srv 25.07.17 }
            endif
          enddo
        enddo                                             !lk 30.06.2017}
#endif
      endif                                               !srv 25.07.17
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
