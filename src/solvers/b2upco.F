*-----------------------------------------------------------------------
*.specification

      subroutine b2upco (nx, ny, nregionv, ua_solve, solvereg,
     &  rxf, pcm0, vol, corpb, pccb0, pb, nb, ub, isb)
      use b2mod_indirect
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_numerics_namelist
     , , only : corr_core_dn
      use b2mod_subsys
      use b2mod_plasma                                                  !srv 07.07.21
     & , only : uap
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, nregionv, isb
      real (kind=R8) ::
     *  rxf, pcm0, vol(-1:nx,-1:ny), corpb(-1:nx,-1:ny),
     *  pccb0(-1:nx,-1:ny)
      logical ua_solve, solvereg(0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *  pb(-1:nx,-1:ny), nb(-1:nx,-1:ny), ub(-1:nx,-1:ny)

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
      integer ncall, ix, iy, ireg
      save ncall
      !$OMP THREADPRIVATE (ncall)
*     ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(-1:nx,-1:ny)
      integer, save :: nrings, jxi, jxa, jsep
      real (kind=R8), save :: na_min = 1.0e4_R8
      !$OMP THREADPRIVATE(na_min, nrings, jxi, jxa, jsep)
      real (kind=R8) :: t0,
     *  pb_(-1:nx,-1:ny), nb_(-1:nx,-1:ny),                               !lk 30.06.17
     *  tt0(-1:ny), tt1(-1:ny), tt2(-1:ny)                                !lk 30.06.17
      logical too_low_density                                             ! IYS 27.11.2017
*   ..procedures
      external xertst, ipgetr, b2xxgs
*   ..initialisation
      data ncall/0/, nrings/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2upco')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call get_jsep(nx,ny,jxi,jxa,jsep)
       call ipgetr ('b2mndr_na_min', na_min)                              ! IYS 27.11.2017
       call ipgeti ('nrings_for_no_speedup_averaging', nrings)            ! IYS 27.03.2019
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test rxf, pcm0
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     .  'faulty argument rxf')
      call xertst (0.0_R8.le.pcm0, 'faulty argument pcm0')

* ..implement the correction
*   ..compute wrk0
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         wrk0(ix,iy) = 0.0_R8
        else
         wrk0(ix,iy) = corpb(ix,iy)-corpb(leftix(ix,iy),leftiy(ix,iy))
        endif
       enddo
      enddo
      call b2xxgs (nx, ny, 0.0_R8, wrk0, 1)
*   ..compute new pb and nb
      do iy = -1, ny
       do ix = -1, nx
        if (0.0_R8.le.corpb(ix,iy)) then
         t0 = 1.0_R8+rxf*corpb(ix,iy)/pb(ix,iy)
        else
         t0 = 1.0_R8/(1.0_R8-rxf*corpb(ix,iy)/pb(ix,iy))
        endif
        nb_(ix,iy) = t0*nb(ix,iy)                        !lk 30.06.2017
        pb_(ix,iy) = t0*pb(ix,iy)                        !lk 30.06.2017
       enddo
      enddo
      if(corr_core_dn(isb).ne.1.0_R8 .and.
     &               (.not. is_neutral(isb))) then               !srv 25.07.17 {  IYS 19.12.2017
*   ..smooth the correction for inner flux surfaces      !lk 30.06.2017{
        tt0 = 0.0e0_R8
        tt1 = 0.0e0_R8
        tt2 = 0.0e0_R8
        do iy = -1,ny
          do ix = -1,nx
            if(on_closed_surface(ix,iy) .and.
     *         iy .le. jsep - nrings .and. bottomiy(ix,iy).ne.-2) then  ! IYS 27.03.2019
              tt0(iy) = tt0(iy) + vol(ix,iy)
              tt1(iy) = tt1(iy) + vol(ix,iy)*(nb_(ix,iy)-nb(ix,iy))
              tt2(iy) = tt2(iy) + vol(ix,iy)*(pb_(ix,iy)-pb(ix,iy))
            endif
          enddo
        enddo
        do iy = -1,ny
          too_low_density = .FALSE.                       ! IYS 27.11.2017
! Check whether the density and/or pressure falls below the threshold 'b2mndr_na_min'
! If it happens, cancel the smoothing of the correction for the whole
! flux surface where such event happens.
          do ix = -1,nx
            if(on_closed_surface(ix,iy) .and.
     *         iy .le. jsep - nrings .and. bottomiy(ix,iy).ne.-2) then    ! IYS 27.03.2019
              t0 = pb(ix,iy)+tt2(iy)/tt0(iy) +
     *         corr_core_dn(isb)*(pb_(ix,iy)-pb(ix,iy)-tt2(iy)/tt0(iy))   ! IYS 19.12.2017
              if (t0 .lt. pb(ix,iy)*na_min / nb(ix,iy)) then
                too_low_density = .TRUE.
                exit
              else
                pb(ix,iy) = t0
              endif
              t0 = nb(ix,iy)+tt1(iy)/tt0(iy) +
     *         corr_core_dn(isb)*(nb_(ix,iy)-nb(ix,iy)-tt1(iy)/tt0(iy))   ! IYS 19.12.2017
              if (t0 .lt. na_min) then
                too_low_density = .TRUE.
                exit
              else
                nb(ix,iy) = t0
              endif
            else                                          !srv 25.07.17 {
              nb(ix,iy) = nb_(ix,iy)
              pb(ix,iy) = pb_(ix,iy)                      !srv 25.07.17 }
            endif
          enddo
          if (too_low_density) then  ! cancel smoothing of the correction on the surface if the corrected density is too low
            do ix=-1,nx
              nb(ix,iy) = nb_(ix,iy)
              pb(ix,iy) = pb_(ix,iy)                      !srv 25.07.17 }
            enddo
          endif
        enddo                                             !lk 30.06.2017}
      else         ! do not smooth the correction
        nb = nb_
        pb = pb_
      endif                                                !srv 25.07.17 }
*   ..compute new ub
      if (ua_solve) then
       if(pcm0*rxf.ne.0.0_R8 .and. all(pccb0.ne.0.0_R8))
     &    uap(:,:,isb) = ub                                             !srv 07.07.21
       do iy = -1, ny
        do ix = -1, nx
*       (note the further damping on the velocity correction)
         ireg = region(ix,iy,0)
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &      solvereg(ireg))
     &    ub(ix,iy) = ub(ix,iy)-pcm0*rxf*pccb0(ix,iy)*vol(ix,iy)*
     &    (wrk0(ix,iy)/(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &     wrk0(rightix(ix,iy),rightiy(ix,iy))/(vol(ix,iy)+
     &      vol(rightix(ix,iy),rightiy(ix,iy))))
        enddo
       enddo
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2upco

      end subroutine b2upco

