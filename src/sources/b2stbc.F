*-----------------------------------------------------------------------
*.specification

      subroutine b2stbc (nCv, nFc, nVx, ns, ismain, ismain0,
     &  switch, geo, mpg, pl, dv, co, rt, st_ext, sr, srw,
     &  psnc, psnl, wrong_flow, main_call)
      use b2mod_types
      use b2mod_constants
#ifndef FIXED_POINT
      use b2mod_balance !djm Jan 2017
     1 , only : b2stbc_sna0to1,b2stbc_smo0to3,b2stbc_she0to3,
     2          b2stbc_shi0to3,balance_netcdf
#endif
      use b2us_feedback
      use b2mod_b2cmpa
      use b2mod_time
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_ad
     & , only : ncall_b2stbc
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_boundary_namelist
     & , only : conpar, enepar, enipar, nbcd, lfeedback, enkpar, potpar,
     &          mompar
      use b2mod_ad
     & , only : my_out_folder,
#ifdef WG_DEBUG
     &          my_out_sf,                                              !srv 27.06.21
#endif
     &          ncall_b2stbc_phys, nsdmax
      use b2mod_math
     & , only : cutlo, cutll, b2mod_math_initialised, small_r4_constant
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, ismain0
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2Source), intent (in) :: sr
      type (B2Coeff), intent (inout) :: co
      type (B2StateExt), intent (in) :: st_ext
      type (B2PlasmaSnapshot), intent (inout) :: psnc, psnl
      logical main_call                                                  !xpb
*   ..output arguments (unspecified on entry)
      type (B2SourceWork), intent (inout) :: srw
      logical wrong_flow                                                 !xpb
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge,
*     due to the boundary conditions.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer i, j, is, iCv, iFc, iFt
      real (kind=R8) :: praverage(2)
      real (kind=R8) :: wrk(nFc)
      real (kind=R8) ::
     &   prsum(2), gssum(2), tesum(2), tisum(2),
     1   nasum(0:ns-1,2),
     1   nef, nif, tef, tif, naf(0:ns-1), zaf(0:ns-1),
     x   nxf, taf, zxf,
     1   volsum(2), qeaverage(2), qiaverage(2), qesum(2), qisum(2)
#ifdef WG_TODO
      integer ixx, regcnt, element_done(-1:nx)
      real (kind=R8) :: t0, testqe, testqi, sum_sch0, sum_vol
#endif
      real (kind=R8) :: intface_s
*   ..procedures
      intrinsic nint
      external xertst, sfill, b2saxpy, ipgetr, ipgeti, intface_s
      external b2xvsg, b2xvps, b2stbc_phys
*   ..initialisation
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc')
*   ..set internal parameters on first call
      if (ncall_b2stbc.eq.0) then
       call xertst (.not.(switch%integral_current.ne.0.0_R8.and.
     .  ((maxval(dv%facdrift).gt.0.0_R8.and.
     .                         switch%b2stbc_fchy_dia.gt.0.0_R8)
     .   .or.switch%b2stbc_neoclassical.gt.0.0_R8)),
     .  'inconsistent boundary settings for current equation')
        lfeedback=switch%b2stbc_feedback.ne.0                          ! dpc-2003-10-23
      endif
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main neutral species ismain; must not be neutral')
      call xertst (0.le.ismain0.and.ismain0.lt.ns,
     . 'invalid main neutral species index ismain0')
      call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     . 'invalid main neutral species ismain0; must be neutral')
*   ..extensive tests on first few calls
      if (ncall_b2stbc.lt.3) then
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
       call b2xvsg (nFc*2, geo%fcVol, 1, 'vol', '.gt.')
*    ..test state
       call b2xvps (nCv, nFc, ns, pl, dv)
*    ..test sign of chce, chci
       wrk(:) = co%chce(:,0)*geo%fcQalf(:,0)
       call b2xvsg (nFc,wrk,1,'chce0','.ge.')
       wrk(:) = co%chce(:,1)*geo%fcQalf(:,1)
       call b2xvsg (nFc,wrk,1,'chce1','.ge.')
       wrk(:) = co%chci(:,0)*geo%fcQalf(:,0)
       call b2xvsg (nFc,wrk,1,'chci0','.ge.')
       wrk(:) = co%chci(:,1)*geo%fcQalf(:,1)
       call b2xvsg (nFc,wrk,1,'chci1','.ge.')
      endif
*   ..test facdrift           !xpb
      call xertst (0.0_R8.le.minval(dv%facdrift) .and.
     & 1.0_R8.ge.maxval(dv%facdrift), 'faulty argument facdrift')
!      if (ncall_b2stbc.eq.1 .and. main_call) then                               !xpb !srv 18.05.09
!        if (maxval(facdrift).eq.0.0_R8) then                             !xpb
!          write (*,*) 'b2stbc_neoclassical =', b2stbc_neoclassical       !xpb
!        else                                                             !xpb
!          write (*,*) 'b2stbc_fchy_dia =', b2stbc_fchy_dia               !xpb
!        endif                                                            !xpb
!      endif
*   ..check kinetic secondary electron emission model
      if (ncall_b2stbc.eq.0 .and. switch%secmodel.eq.1) then
        call xertst (nint(zn(ismain)).eq.1,
     &   'Kinetic secondary emission model only for hydrogen plasmas!')
        write (*,*) 'Kinetic secondary emission model turned on'
        write (*,*) 'Graphite target plates are assumed'
        if (nint(zamax(ismain)).ne.1) write(*,*)
     .   'WARNING! Kinetic secondary emission expects HYDROGEN plasmas!'
      endif
*   ..initialise sources to 0
      call sfill (nCv*2*ns, 0.0_R8, srw%sna0, 1)
      call sfill (nCv*4*ns, 0.0_R8, srw%smo0, 1)
      call sfill (nCv*2, 0.0_R8, srw%sne0, 1)
      call sfill (nCv*4, 0.0_R8, srw%she0, 1)
      call sfill (nCv*4, 0.0_R8, srw%shi0, 1)
      call sfill (nCv*4, 0.0_R8, srw%shn0, 1)
      call sfill (nCv*4, 0.0_R8, srw%sch0, 1)
      call sfill (nCv*4, 0.0_R8, srw%skt0, 1)
      call sfill (nCv*4, 0.0_R8, srw%szt0, 1)

* ..compute standard form volume sources
*   (This code placed here for want of a better location)
*   ..low-level particle source
      do is = 0, ns-1
       call b2saxpy (nCv,switch%sna0ep,geo%cvVol,1,srw%sna0(1,0,is),1)
      enddo
      call b2saxpy (nCv, switch%she0ep, geo%cvVol, 1, srw%she0(1,0), 1)
      call b2saxpy (nCv, switch%shi0ep, geo%cvVol, 1, srw%shi0(1,0), 1)

*   ..find the average pressure for inner flux surface           !sv 24.03.99
      iFt = mpg%cvFt(omp(1))                                     !sv 24.03.99
      tesum(1)=0.0_R8
      tisum(1)=0.0_R8
      nasum(:,1)=0.0_R8
      prsum(1)=0.0_R8                                                    !xpb
      qesum(1)=0.0_R8
      qisum(1)=0.0_R8
      gssum(1)=0.0_R8                                                    !xpb
      volsum(1)=0.0_R8
      do i = 1, mpg%ftCvP(iFt,2)                                         !xpb
        iCv = mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
        if(mpg%cvOnClosedSurface(iCv)) then                              !xpb
          do j = 1, mpg%cvFcP(iCv,2)
            iFc = mpg%cvFc(mpg%cvFcP(iCv,1)+j-1)
            if ((mpg%fcCv(iFc,1).eq.iCv.and.mpg%fcCv(iFc,2).gt.mpg%nCi)
     &          .or.
     &          (mpg%fcCv(iFc,2).eq.iCv.and.mpg%fcCv(iFc,1).gt.mpg%nCi))
     >       then
              tef = intface_s(iFc,nCv,nFc,mpg%fcCv,geo%fcVol,pl%te)
              tif = intface_s(iFc,nCv,nFc,mpg%fcCv,geo%fcVol,pl%ti)
              nef = intface_s(iFc,nCv,nFc,mpg%fcCv,geo%fcVol,dv%ne)
              nif = intface_s(iFc,nCv,nFc,mpg%fcCv,geo%fcVol,dv%ni(1,1))
              do is = 0, ns-1
                naf(is) =
     &           intface_s(iFc,nCv,nFc,mpg%fcCv,geo%fcVol,pl%na(1,is))
                zaf(is) =
     &           intface_s(iFc,nCv,nFc,mpg%fcCv,geo%fcVol,rt%rza(1,is))
              end do
              tesum(1)=tesum(1) + tef*geo%fcS(iFc)
              tisum(1)=tisum(1) + tif*geo%fcS(iFc)
              nasum(:,1)=nasum(:,1) + naf(:)*geo%fcS(iFc)
              prsum(1)=prsum(1) +
     1         (nef*tef+nif*tif)*geo%cvVol(iCv)*
     2          geo%fcBb(iFc,2)*geo%fcHz(iFc)
              qesum(1)=qesum(1)+2.5_R8*geo%fcHz(iFc)*
     1          nef*tef**2*geo%fcBb(iFc,2)/qe*geo%fcS(iFc)
              do is = 0, ns-1
                if (.not.is_neutral(is)) then
                  qisum(1)=qisum(1)-2.5_R8*geo%fcHz(iFc)*
     3              naf(is)*tif**2*geo%fcBb(iFc,2)/
     5             (qe*zaf(is))*geo%fcS(iFc)
                endif
              enddo
              do is = 0, st_ext%ns-1
                nxf = intface_s(iFc,nCv,nFc,mpg%fcCv,geo%fcVol,
     &           st_ext%na(1,is))
                taf = intface_s(iFc,nCv,nFc,mpg%fcCv,geo%fcVol,
     &           st_ext%ta(1,is))
                prsum(1)=prsum(1) + nxf*taf*geo%cvVol(iCv)*
     2           geo%fcBb(iFc,2)*geo%cvHz(iCv)
                if (.not.st_ext%is_neutral(is)) then
                  zxf = intface_s(iFc,nCv,nFc,mpg%fcCv,geo%fcVol,
     &             st_ext%za(1,is))
                  qisum(1)=qisum(1)-2.5_R8*geo%cvHz(iFc)*
     3             nxf*taf**2*geo%fcBb(iFc,2)/(qe*zxf)*geo%fcS(iFc)
                endif
              enddo
              gssum(1)=gssum(1) + geo%fcS(iFc)                           !xpb
              volsum(1)=volsum(1) + geo%cvVol(iCv)
            endif                                                        !xpb
          enddo
        endif
      enddo                                                      !sv 24.03.99
      if (volsum(1).gt.0.0_R8) praverage(1)=prsum(1)/volsum(1)           !xpb
      if (gssum(1).gt.0.0_R8) then
        qeaverage(1)=qesum(1)/gssum(1)
        qiaverage(1)=qisum(1)/gssum(1)
      end if
#ifdef WG_TODO
*   ..find the average pressure for island center flux surface           !xpb
      iy = ny                                                            !xpb
      tesum(2)=0.0_R8
      tisum(2)=0.0_R8
      nasum(:,2)=0.0_R8
      prsum(2)=0.0_R8                                                    !xpb
      qesum(2)=0.0_R8
      qisum(2)=0.0_R8
      gssum(2)=0.0_R8                                                    !xpb
      volsum(2)=0.0_R8
      praverage(2)=0.0_R8                                                !xpb
      qeaverage(2)=0.0_R8                                                !xpb
      qiaverage(2)=0.0_R8                                                !xpb
      do ix = -1, nx                                                     !xpb
        if(on_closed_surface(ix,iy)) then                                !xpb
          tef = (te(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           te(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          tif = (ti(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           ti(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          nef = (ne(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           ne(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          nif = (ni(ix,iy,1)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           ni(bottomix(ix,iy),bottomiy(ix,iy),1)*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          naf(:) = (na(ix,iy,:)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           na(bottomix(ix,iy),bottomiy(ix,iy),:)*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          zaf(:) = (rza(ix,iy,:)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           rza(bottomix(ix,iy),bottomiy(ix,iy),:)*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          tesum(2)=tesum(2) + tef*gs(ix,iy,1)                            !xpb
          tisum(2)=tisum(2) + tif*gs(ix,iy,1)                            !xpb
          nasum(:,2)=nasum(:,2) + naf(:)*gs(ix,iy,1)                     !xpb
          prsum(2)=prsum(2) +
     1     (nef*tef+nif*tif)*vol(ix,iy)*
     2     wbbv(bottomix(ix,iy),bottomiy(ix,iy),2)*
     3     vol(ix,iy)/gs(ix,iy,2)
          qesum(2)=qesum(2)+2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     1     nef*tef**2*wbbv(bottomix(ix,iy),bottomiy(ix,iy),2)/qe*
     2     gs(ix,iy,1)
          do is = 0, ns-1
            if (.not.is_neutral(is)) then
              qisum(2)=qisum(2)-2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     3         naf(is)*tif**2*wbbv(bottomix(ix,iy),bottomiy(ix,iy),2)/
     5         (qe*zaf(is))*gs(ix,iy,1)
            endif
          enddo
          do is = 0, ns_ext-1
            nxf = (na_ext(ix,iy,is)*
     .             hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &             na_ext(bottomix(ix,iy),bottomiy(ix,iy),is)*
     .             hy1(ix,iy))/
     &            (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
            taf = (ta_ext(ix,iy,is)*
     .             hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &             ta_ext(bottomix(ix,iy),bottomiy(ix,iy),is)*
     .             hy1(ix,iy))/
     &            (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
            prsum(2)=prsum(2) + nxf*taf*vol(ix,iy)
     2        *wbbv(bottomix(ix,iy),bottomiy(ix,iy),2)
     3        *vol(ix,iy)/gs(ix,iy,2)
            if (.not.is_ext_neutral(is)) then
              zxf = (za_ext(ix,iy,is)*
     .               hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &               za_ext(bottomix(ix,iy),bottomiy(ix,iy),is)*
     .               hy1(ix,iy))/
     &              (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
              qisum(2)=qisum(2)-2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     3         nxf*taf**2*wbbv(bottomix(ix,iy),bottomiy(ix,iy),2)/
     5         (qe*zxf)*gs(ix,iy,1)
            endif
          enddo
          gssum(2)=gssum(2) + gs(ix,iy,1)                                !xpb
          volsum(2)=volsum(2) + vol(ix,iy)
        endif                                                            !xpb
      enddo                                                              !xpb
      if (volsum(2).gt.0.0_R8) praverage(2)=prsum(2)/volsum(2)           !xpb
      if (gssum(2).gt.0.0_R8) then
        qeaverage(2)=qesum(2)/gssum(2)                                   !xpb
        qiaverage(2)=qisum(2)/gssum(2)                                   !xpb
      endif
#endif
      if(switch%b2stbc_boundary_namelist.ge.1) then
        if(lfeedback.or.switch%b2stbc_feedback.ne.0)
     &    call compute_feedback(nCv, nFc, ns, ismain,
     &     switch, geo, mpg, pl, dv, rt, sr, srw, psnc, psnl,
     &     main_call)

        if(ncall_b2stbc.eq.0) write(*,*) ' b2stbc_phys called'
        call b2stbc_phys(nCv, nFc, nVx, ns, ismain, ismain0,
     &   switch, geo, mpg, pl, dv, co, rt, st_ext, srw,
     &   wrong_flow, main_call, praverage)
        if(wrong_flow) write(*,*)
     &        'b2stbc: wrong_flow returned from b2stbc_phys'
      endif

!WG_RM      if(lfeedback.or.switch%b2stbc_feedback.ne.0)
!WG_RM     &  call b2stbc_fb_post(nCv, nFc, ns, ismain,
!WG_RM     &    switch, geo, pl, dv, srw, main_call)

#ifdef WG_TODO
c added 2001.11.26 dpc/xpb
c to fix problems with the core energy boundary condition arising
c from diamagnetic drifts
      testqe=0.0_R8
      testqi=0.0_R8
      if(maxval(facdrift).ne.0.0_R8) then
        do ix=-1,nx
          iy=0                                        ! south side
          if(.not.on_closed_surface(ix,iy) .and.
     &       b2stbc_fchy_dia_coreonly.eq.0) then  ! code for private flux
            do is=0,ns-1
              if(.not.is_neutral(is)) then
                t0 =
     &           -2.5_R8*
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
     &           /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &           *(rza(bottomix(ix,iy),bottomiy(ix,iy),is)
     &            *na(bottomix(ix,iy),bottomiy(ix,iy),is)
     &            *te(bottomix(ix,iy),bottomiy(ix,iy))
     &            /ti(bottomix(ix,iy),bottomiy(ix,iy))
     &           +rza(ix,iy,is)*na(ix,iy,is)*te(ix,iy)/ti(ix,iy))/
     1           2.0_R8*vadia(ix,iy,1,1,is)*b2stbc_fchy_dia*
     1           (te(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &            te(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &           (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
                she0(ix,-1,0) = she0(ix,-1,0) + t0
                testqe = testqe + t0
                t0 =
     &           +2.5_R8*
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
     &           /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &           *(na(bottomix(ix,iy),bottomiy(ix,iy),is)
     &            +na(ix,iy,is))/2.0_R8
     &           *vadia(ix,iy,1,1,is)*b2stbc_fchy_dia*
     1           (ti(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &            ti(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &           (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
                shi0(ix,-1,0) = shi0(ix,-1,0) + t0
                testqi = testqi + t0
              endif
            enddo
          elseif(on_closed_surface(ix,iy)) then
            if (diagno.ge.2) write(*,*) 'Before she0,shi0 =',
     &            ix,she0(ix,-1,0),shi0(ix,-1,0)
            she0(ix,-1,0)=she0(ix,-1,0)+qeaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            testqe=testqe+qeaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            shi0(ix,-1,0)=shi0(ix,-1,0)+qiaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            testqi=testqi+qiaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            if (diagno.ge.2) write(*,*) 'Electron term =',
     .       +qeaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            if (diagno.ge.2) write(*,*) 'Ion term =',
     .       +qiaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
          endif
          iy = ny
          if(on_closed_surface(ix,iy)) then  ! code for island center
            she0(ix,iy,0)=she0(ix,iy,0)-qeaverage(2)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            shi0(ix,iy,0)=shi0(ix,iy,0)-qiaverage(2)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
          elseif(region(ix,iy,0).ne.5 .and.
     &           b2stbc_fchy_dia_coreonly.eq.0) then  ! North boundary
            do is = 0, ns-1
              if (.not.is_neutral(is)) then
                t0 =
     &           +2.5_R8*
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
     &           /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &           *(rza(bottomix(ix,iy),bottomiy(ix,iy),is)
     &           *na(bottomix(ix,iy),bottomiy(ix,iy),is)
     &           *te(bottomix(ix,iy),bottomiy(ix,iy))
     &           /ti(bottomix(ix,iy),bottomiy(ix,iy))
     &           +rza(ix,iy,is)*na(ix,iy,is)*te(ix,iy)/ti(ix,iy))/
     1           2.0_R8*vadia(ix,iy,1,1,is)*b2stbc_fchy_dia*
     1           (te(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &            te(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &           (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
                she0(ix,iy,0)=she0(ix,iy,0)+t0
                testqe=testqe+t0
                t0 =
     &           -2.5_R8*
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
     &           /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &           *(na(bottomix(ix,iy),bottomiy(ix,iy),is)
     &            +na(ix,iy,is))/2.0_R8
     &           *vadia(ix,iy,1,1,is)*b2stbc_fchy_dia*
     1           (ti(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &            ti(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &           (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
                shi0(ix,iy,0)=shi0(ix,iy,0)+t0
                testqi=testqi+t0
              endif
            enddo
          endif
        enddo
      endif
      if (diagno.ge.2) write(*,*) 'testqe, testqi ', testqe, testqi
c end 2001.11.26 dpc/xpb addition

c corrections for N & S current boundary conditions
c if we want the integral current to be 0
      if(integral_current.gt.0.0_R8) then
        regcnt=0
        do iy=-1,ny,ny+1                         ! do S + N boundaries
          element_done=0
          ix=-1
          do while(ix.le.nx)                     ! do each ix
            if(element_done(ix).eq.0) then       ! have we done this chain yet?
              regcnt=regcnt+1
              sum_sch0=sch0(ix,iy,0)+
     +         sch0(ix,iy,1)*po(ix,iy)+
     +         sch0(ix,iy,2)*ne(ix,iy)+
     +         sch0(ix,iy,3)*ne(ix,iy)*po(ix,iy)
              sum_vol=vol(ix,iy)
              element_done(ix)=regcnt
              ixx=rightix(ix,iy)
              do while(ixx.ge.-1.and.ixx.le.nx)  ! do chain of elements connected to this element
                if(element_done(ixx).ne.0) then
                  ixx=nx+1                       ! get here for periodic case
                else
                  sum_sch0=sum_sch0+sch0(ixx,iy,0)+
     +             sch0(ixx,iy,1)*po(ixx,iy)+
     +             sch0(ixx,iy,2)*ne(ixx,iy)+
     +             sch0(ixx,iy,3)*ne(ixx,iy)*po(ixx,iy)
                  sum_vol=sum_vol+vol(ixx,iy)
                  element_done(ixx)=regcnt
                  ixx=rightix(ixx,iy)
                endif
              enddo
              write(*,*) 'Subtracting ',sum_sch0,' for region ',regcnt
              do ixx=-1,nx
                if(element_done(ixx).eq.regcnt)
     1           sch0(ixx,iy,0)=sch0(ixx,iy,0)-
     2           sum_sch0*vol(ixx,iy)/sum_vol*integral_current
              enddo
            endif
            ix=ix+1
          enddo
          if (ncall_b2stbc.eq.0 .or. diagno.ge.2)
     >     write(*,'(i5/(10(1x,10i1)))') iy,element_done
        enddo
      endif
#endif
* ..return
#ifdef WG_TODO
      if((ncall_b2stbc.eq.0 .and. diagno.ge.2) .or. diagno.ge.3) then
        do ix=-1,nx
          write(*,'(a,i4,1p,6(1x,e14.7))') 'S:ix,ti,ni,shi0 ',
     1     ix,ti(ix,-1)/ev,ni(ix,-1,0),(shi0(ix,-1,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'S:ix,te,ne,she0 ',
     1     ix,te(ix,-1)/ev,ne(ix,-1),(she0(ix,-1,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'N:ix,ti,ni,shi0 ',
     1     ix,ti(ix,ny)/ev,ni(ix,ny,0),(shi0(ix,ny,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'N:ix,te,ne,she0 ',
     1     ix,te(ix,ny)/ev,ne(ix,ny),(she0(ix,ny,is),is=0,3)
        enddo
        do iy=-1,ny
          write(*,'(a,i4,1p,6(1x,e14.7))') 'W:iy,ti,ni,shi0 ',
     1     iy,ti(-1,iy)/ev,ni(-1,iy,0),(shi0(-1,iy,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'W:iy,te,ne,she0 ',
     1     iy,te(-1,iy)/ev,ne(-1,iy),(she0(-1,iy,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'E:iy,ti,ni,shi0 ',
     1     iy,ti(nx,iy)/ev,ni(nx,iy,0),(shi0(nx,iy,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'E:iy,te,ne,she0 ',
     1     iy,te(nx,iy)/ev,ne(nx,iy),(she0(nx,iy,is),is=0,3)
        enddo
      endif
#endif
c      write(*,*) 'Copying sna0 to b2stbc_sna'
      srw%b2stbc_sna(:,0:ns-1)=
     1 srw%sna0(:,0,0:ns-1)+
     2 srw%sna0(:,1,0:ns-1)*pl%na(:,0:ns-1)
c      write(*,*) 'Copying smo0 to b2stbc_smo'
      do is=0,ns-1
        srw%b2stbc_smo(:,is)=
     1   srw%smo0(:,0,is)+
     2   srw%smo0(:,1,is)*pl%ua(:,is)+
     3   srw%smo0(:,2,is)*am(is)*mp*pl%na(:,is)+
     4   srw%smo0(:,3,is)*am(is)*mp*pl%na(:,is)*pl%ua(:,is)
      enddo
c      write(*,*) 'Copying she0 to b2stbc_she'
      srw%b2stbc_she(:)=
     1 srw%she0(:,0)+
     2 srw%she0(:,1)*pl%te(:)+
     3 srw%she0(:,2)*dv%ne(:)+
     4 srw%she0(:,3)*dv%ne(:)*pl%te(:)
c      write(*,*) 'Copying shi0 to b2stbc_shi'
      srw%b2stbc_shi(:)=
     1 srw%shi0(:,0)+
     2 srw%shi0(:,1)*pl%ti(:)+
     3 srw%shi0(:,2)*dv%ni(:,0)+
     4 srw%shi0(:,3)*dv%ni(:,0)*pl%ti(:)
      srw%b2stbc_shn(:)=
     1 srw%shn0(:,0)+
     2 srw%shn0(:,1)*pl%tn(:)+
     3 srw%shn0(:,2)*dv%nn(:)+
     4 srw%shn0(:,3)*dv%nn(:)*pl%tn(:)
c      write(*,*) 'Copying sch0 to b2stbc_sch'
      srw%b2stbc_sch(:)=
     1 srw%sch0(:,0)+
     2 srw%sch0(:,1)*pl%po(:)+
     3 srw%sch0(:,2)*dv%ne(:)+
     4 srw%sch0(:,3)*dv%ne(:)*pl%po(:)
      srw%b2stbc_skt(:)=
     1 srw%skt0(:,0)+
     2 srw%skt0(:,1)*pl%kt(:)+
     3 srw%skt0(:,2)*dv%ni(:,1)+
     4 srw%skt0(:,3)*dv%ni(:,1)*pl%kt(:)
      srw%b2stbc_szt(:)=
     1 srw%szt0(:,0)+
     2 srw%szt0(:,1)*pl%zt(:)+
     3 srw%szt0(:,2)*dv%ni(:,1)+
     4 srw%szt0(:,3)*dv%ni(:,1)*pl%zt(:)
c      write(*,*) 'Copying sne0 to b2stbc_sne'
      srw%b2stbc_sne(:)=
     1 srw%sne0(:,0)+
     2 srw%sne0(:,1)*dv%ne(:)
c      write(*,*) 'Copied s..0 to b2stbc_s..'

#ifndef FIXED_POINT
cdjm Jan2017 Keep linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2stbc_sna0to1 = srw%sna0
        b2stbc_smo0to3 = srw%smo0
        b2stbc_she0to3 = srw%she0
        b2stbc_shi0to3 = srw%shi0
      endif
#endif

      if(ncall_b2stbc.eq.0 .or. main_call) ncall_b2stbc = ncall_b2stbc+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbc

      end subroutine b2stbc
