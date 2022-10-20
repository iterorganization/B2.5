*-----------------------------------------------------------------------
*.specification

      subroutine b2stbc (nCv, nFc, nVx, ns, ismain, ismain0, BoRiS,
     &  switch, geo, mpg, pl, dv, co, rt, rtw, st_ext, srw,
     &  psnc, psnl, wrong_flow, main_call)
      use b2mod_types
!      use b2mod_geo_corner
      use b2mod_constants
!      use b2mod_sources
!     1 , only : b2stbc_sna,b2stbc_smo,b2stbc_she,b2stbc_shi,b2stbc_sch,
!     2          b2stbc_sne
!      use b2mod_balance !djm Jan 2017
!     1 , only : b2stbc_sna0to1,b2stbc_smo0to3,b2stbc_she0to3,
!     2          b2stbc_shi0to3,balance_netcdf
      use b2us_feedback
      use b2mod_b2cmpa
      use b2mod_time
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
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
     &          ncall_b2stbc,
     &          ncall_b2stbc_phys, nsdmax
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, ismain0
      real (kind=R8) ::
     *  BoRiS
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2RatesWork), intent (in) :: rtw
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
      integer iCv, is, regcnt, ixx
      real (kind=R8) ::
     *   t0, testqe, testqi,
     &   praverage(2), prsum(2), gssum(2), tesum(2), tisum(2),
     1   nasum(0:ns-1,2),
     1   nef, nif, tef, tif, naf(0:ns-1), zaf(0:ns-1),
     x   nxf, taf, zxf,
     1   volsum(2), qeaverage(2), qiaverage(2), qesum(2), qisum(2),
     1   sum_sch0, sum_vol
!WG_TODO      integer element_done(-1:nx)
*   ..procedures
      intrinsic nint
      real (kind=R8) ::
     *   b2ssum
      external xertst, sfill, b2ssum, b2saxpy, ipgetr, ipgeti
      external b2xvsg, b2xvff, b2xvps, b2xxgs,
     1   b2stbc_phys
#ifdef WG_TODO
     &          b2stbc_bas                                              !srv 27.06.21
#endif
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
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
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
       call b2xvsg (nCv,co%chce(:,0)*geo%fcQalf(:,0),1,'chce0','.ge.')
       call b2xvsg (nCv,co%chce(:,1)*geo%fcQalf(:,1),1,'chce1','.ge.')
       call b2xvsg (nCv,co%chci(:,0)*geo%fcQalf(:,0),1,'chci0','.ge.')
       call b2xvsg (nCv,co%chci(:,1)*geo%fcQalf(:,1),1,'chci1','.ge.')
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

#ifdef WG_TODO
*   ..find the average pressure for inner flux surface           !sv 24.03.99
      iy = -1                                                    !sv 24.03.99
      tesum(1)=0.0_R8
      tisum(1)=0.0_R8
      nasum(:,1)=0.0_R8
      prsum(1)=0.0_R8                                                    !xpb
      qesum(1)=0.0_R8
      qisum(1)=0.0_R8
      gssum(1)=0.0_R8                                                    !xpb
      volsum(1)=0.0_R8
      do ix = -1, nx                                                     !xpb
        if(on_closed_surface(ix,iy)) then                                !xpb
          tef = (te(ix,iy)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           te(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          tif = (ti(ix,iy)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           ti(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          nef = (ne(ix,iy)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           ne(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          nif = (ni(ix,iy,1)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           ni(topix(ix,iy),topiy(ix,iy),1)*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          naf(:) = (na(ix,iy,:)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           na(topix(ix,iy),topiy(ix,iy),:)*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          zaf(:) = (rza(ix,iy,:)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           rza(topix(ix,iy),topiy(ix,iy),:)*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          tesum(1)=tesum(1) + tef*gs(topix(ix,iy),topiy(ix,iy),1)
          tisum(1)=tisum(1) + tif*gs(topix(ix,iy),topiy(ix,iy),1)
          nasum(:,1)=nasum(:,1) + naf(:)*gs(topix(ix,iy),topiy(ix,iy),1)
          prsum(1)=prsum(1) +
     1     (nef*tef+nif*tif)*vol(ix,iy)
     2     *wbbv(topix(ix,iy),topiy(ix,iy),2)*vol(ix,iy)/gs(ix,iy,2)
          qesum(1)=qesum(1)+2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     1     nef*tef**2*wbbv(topix(ix,iy),topiy(ix,iy),2)/qe*
     2     gs(topix(ix,iy),topiy(ix,iy),1)
          do is = 0, ns-1
            if (.not.is_neutral(is)) then
              qisum(1)=qisum(1)-2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     3         naf(is)*tif**2*wbbv(topix(ix,iy),topiy(ix,iy),2)/
     5         (qe*zaf(is))*gs(topix(ix,iy),topiy(ix,iy),1)
            endif
          enddo
          do is = 0, ns_ext-1
            nxf = (na_ext(ix,iy,is)*hy1(topix(ix,iy),topiy(ix,iy))+
     &             na_ext(topix(ix,iy),topiy(ix,iy),is)*hy1(ix,iy))/
     &            (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
            taf = (ta_ext(ix,iy,is)*hy1(topix(ix,iy),topiy(ix,iy))+
     &             ta_ext(topix(ix,iy),topiy(ix,iy),is)*hy1(ix,iy))/
     &            (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
            prsum(1)=prsum(1) + nxf*taf*vol(ix,iy)
     2        *wbbv(topix(ix,iy),topiy(ix,iy),2)*vol(ix,iy)/gs(ix,iy,2)
            if (.not.is_ext_neutral(is)) then
              zxf = (za_ext(ix,iy,is)*hy1(topix(ix,iy),topiy(ix,iy))+
     &               za_ext(topix(ix,iy),topiy(ix,iy),is)*hy1(ix,iy))/
     &              (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
              qisum(1)=qisum(1)-2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     3         nxf*taf**2*wbbv(topix(ix,iy),topiy(ix,iy),2)/
     5         (qe*zxf)*gs(topix(ix,iy),topiy(ix,iy),1)
            endif
          enddo
          gssum(1)=gssum(1) + gs(topix(ix,iy),topiy(ix,iy),1)            !xpb
          volsum(1)=volsum(1) + vol(ix,iy)
        endif                                                            !xpb
      enddo                                                      !sv 24.03.99
      if (volsum(1).gt.0.0_R8) praverage(1)=prsum(1)/volsum(1)           !xpb
      if (gssum(1).gt.0.0_R8) then
        qeaverage(1)=qesum(1)/gssum(1)
        qiaverage(1)=qisum(1)/gssum(1)
      end if
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
      if (ncall_b2stbc.eq.0) lbndusr=.false.                                    !xpb
#endif
      if(switch%b2stbc_boundary_namelist.ge.1) then
        if(lfeedback.or.switch%b2stbc_feedback.ne.0)
     &    call compute_feedback(nCv, nFc, nVx, ns, ismain,
     &     switch, geo, mpg, pl, dv, rt, srw, psnc, psnl, 
     &     main_call, switch%b2stbc_diagno)

        if(ncall_b2stbc.eq.0) write(*,*) ' b2stbc_phys called'
        call b2stbc_phys(nCv, nFc, nVx, ns, ismain, ismain0, BoRiS,
     &   switch, geo, mpg, pl, dv, co, rt, rtw, st_ext, srw,
     &   wrong_flow, main_call, switch%b2stbc_diagno,
     &   switch%b2stbc_fchy_dia, switch%b2stbc_neoclassical, praverage,
     &   switch%b2stbc_fchy_dia_coreonly,
     &   switch%secmodel, switch%sheath_drift_fix)
        if(wrong_flow) write(*,*)
     &        'b2stbc: wrong_flow returned from b2stbc_phys'
#ifdef WG_TODO
        if(lbndusr) then
          if(ncall_b2stbc.eq.0) write(*,*) 'BNDUSR called'
          call b2stbc_bas(nx, ny, ns, ismain,
     1     facdrift, vol, hx, hy, qz, qc, gs,
     1     pbs, bb, na, ua, uadia, vaecrb, te, ti, po, ne, ni,
     1     fna, fhe, fhi, fch, fne, fni, chce, chve, chci, chvi, calf,
     &     sna0, smo0, she0, shi0, sch0, sne0,
     &     wrong_flow, main_call, diagno,
     &     b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     &     b2stbc_fchy_dia_coreonly, secmodel, lbndusr,
     &     area_fix, sheath_drift_fix)
          if(wrong_flow) write(*,*)
     &           'b2stbc: wrong_flow returned from b2stbc_bas'
        endif      !lbndusr
#endif
      else      !boundary_namelist.ge.1
#ifdef WG_TODO
        call b2stbc_bas(nx, ny, ns, ismain,
     1   facdrift, vol, hx, hy, qz, qc, gs,
     1   pbs, bb, na, ua, uadia, vaecrb, te, ti, po, ne, ni,
     1   fna, fhe, fhi, fch, fne, fni, chce, chve, chci, chvi, calf,
     &   sna0, smo0, she0, shi0, sch0, sne0,
     &   wrong_flow, main_call, diagno,
     &   b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     &   b2stbc_fchy_dia_coreonly, secmodel, lbndusr,
     &   area_fix, sheath_drift_fix)
        if(wrong_flow) write(*,*)
     &       'b2stbc: wrong_flow returned from b2stbc_bas'
#endif
      endif

!WG_RM      if(lfeedback.or.switch%b2stbc_feedback.ne.0)
!WG_RM     &  call b2stbc_fb_post(nCv, nFc, nVx, ns, ismain, BoRiS,
!WG_RM     &    switch, geo, mpg, pl, dv, srw, main_call, 
!WG_RM     &    switch%b2stbc_diagno)

* ..compute conditions on the cuts
*     (no flow)
!WG_TODO*   ..compute suitable large coefficient
!WG_TODO      t0 = 1.0e20_R8*b2ssum(nCv,geo%cvVol,1)
!WG_TODO*     (the dimension of t0 is m**3/s)
!WG_TODO*   ..loop over species
!WG_TODO      do is = 0, ns-1
!WG_TODO*    ..set parallel momentum source
!WG_TODO       call b2xxgs (nx, ny, 0.0_R8, smo0(-1,-1,0,is), 1)
!WG_TODO       call b2xxgs (nx, ny, 0.0_R8, smo0(-1,-1,1,is), 1)
!WG_TODO       call b2xxgs (nx, ny, 0.0_R8, smo0(-1,-1,2,is), 1)
!WG_TODO       call b2xxgs (nx, ny, -t0, smo0(-1,-1,3,is), 1)
!WG_TODO      enddo

c added 2001.11.26 dpc/xpb
c to fix problems with the core energy boundary condition arising
c from diamagnetic drifts
      testqe=0.0_R8
      testqi=0.0_R8
#ifdef WG_TODO
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
c      write(*,*) 'Copying sna0 to b2stbc_sna'
      b2stbc_sna(-1:nx,-1:ny,0:ns-1)=
     1 sna0(-1:nx,-1:ny,0,0:ns-1)+
     2 sna0(-1:nx,-1:ny,1,0:ns-1)*na(-1:nx,-1:ny,0:ns-1)
c      write(*,*) 'Copying smo0 to b2stbc_smo'
      do is=0,ns-1
        b2stbc_smo(-1:nx,-1:ny,is)=
     1   smo0(-1:nx,-1:ny,0,is)+
     2   smo0(-1:nx,-1:ny,1,is)*ua(-1:nx,-1:ny,is)+
     3   smo0(-1:nx,-1:ny,2,is)*am(is)*mp*na(-1:nx,-1:ny,is)+
     4   smo0(-1:nx,-1:ny,3,is)*am(is)*mp*na(-1:nx,-1:ny,is)*
     5   ua(-1:nx,-1:ny,is)
      enddo
c      write(*,*) 'Copying she0 to b2stbc_she'
      b2stbc_she(-1:nx,-1:ny)=
     1 she0(-1:nx,-1:ny,0)+
     2 she0(-1:nx,-1:ny,1)*te(-1:nx,-1:ny)+
     3 she0(-1:nx,-1:ny,2)*ne(-1:nx,-1:ny)+
     4 she0(-1:nx,-1:ny,3)*ne(-1:nx,-1:ny)*te(-1:nx,-1:ny)
c      write(*,*) 'Copying shi0 to b2stbc_shi'
      b2stbc_shi(-1:nx,-1:ny)=
     1 shi0(-1:nx,-1:ny,0)+
     2 shi0(-1:nx,-1:ny,1)*ti(-1:nx,-1:ny)+
     3 shi0(-1:nx,-1:ny,2)*ni(-1:nx,-1:ny,0)+
     4 shi0(-1:nx,-1:ny,3)*ni(-1:nx,-1:ny,0)*ti(-1:nx,-1:ny)
c      write(*,*) 'Copying sch0 to b2stbc_sch'
      b2stbc_sch(-1:nx,-1:ny)=
     1 sch0(-1:nx,-1:ny,0)+
     2 sch0(-1:nx,-1:ny,1)*po(-1:nx,-1:ny)+
     3 sch0(-1:nx,-1:ny,2)*ne(-1:nx,-1:ny)+
     4 sch0(-1:nx,-1:ny,3)*ne(-1:nx,-1:ny)*po(-1:nx,-1:ny)
c      write(*,*) 'Copying sne0 to b2stbc_sne'
      b2stbc_sne(-1:nx,-1:ny)=
     1 sne0(-1:nx,-1:ny,0)+
     2 sne0(-1:nx,-1:ny,1)*ne(-1:nx,-1:ny)
c      write(*,*) 'Copied s..0 to b2stbc_s..'

cdjm Jan2017 Keep linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2stbc_sna0to1 = sna0
        b2stbc_smo0to3 = smo0
        b2stbc_she0to3 = she0
        b2stbc_shi0to3 = shi0
      endif
#endif
      if(ncall_b2stbc.eq.0 .or. main_call) ncall_b2stbc = ncall_b2stbc+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbc

      end subroutine b2stbc

*-----------------------------------------------------------------------
*.specification
      function bcint(fun,nx,ny,il,ib)
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      use b2mod_boundary_namelist
      implicit none
*   ..output result (undefined upon entry)
      real (kind=R8) :: bcint
*   ..input arguments (unchanged on exit)
      integer nx, ny, il, ib
      real (kind=R8) :: fun(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     BCINT provides interpolated face-centred values for cell-centred
*     quantities (array fun), for use in boundary conditions.
*     The interpolation is done according to distance away from the
*     cell face. The boundary index is ib, and the cell index along
*     the boundary is il.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy
*   ..procedures
      real (kind=R8) ::
     *  hy1
      external xertst, xerrab
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*  ..test input parameters
      call xertst (0.le.nx.and.0.le.ny, 'faulty parameter nx, ny')
      call xertst (0.lt.ib.and.nbc.ge.ib, 'faulty parameter ib')
      call xertst (1.le.il.and.bc_list_size(ib).ge.il,
     .  'faulty parameter il')

* ..main computation
      if (bcchar(ib).eq.'S') then
        bcint =
     .   (fun(bc_list_x(il,ib),bc_list_y(il,ib))*
     .    hy1(topix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        topiy(bc_list_x(il,ib),bc_list_y(il,ib)))+
     .    fun(topix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        topiy(bc_list_x(il,ib),bc_list_y(il,ib)))*
     .    hy1(bc_list_x(il,ib),bc_list_y(il,ib)))
     .  /(hy1(bc_list_x(il,ib),bc_list_y(il,ib))+
     .    hy1(topix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        topiy(bc_list_x(il,ib),bc_list_y(il,ib))))
      elseif (bcchar(ib).eq.'N') then
        bcint =
     .   (fun(bc_list_x(il,ib),bc_list_y(il,ib))*
     .    hy1(bottomix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        bottomiy(bc_list_x(il,ib),bc_list_y(il,ib)))+
     .    fun(bottomix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        bottomiy(bc_list_x(il,ib),bc_list_y(il,ib)))*
     .    hy1(bc_list_x(il,ib),bc_list_y(il,ib)))
     .  /(hy1(bc_list_x(il,ib),bc_list_y(il,ib))+
     .    hy1(bottomix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        bottomiy(bc_list_x(il,ib),bc_list_y(il,ib))))
      elseif (bcchar(ib).eq.'E') then
        bcint =
     .   (fun(bc_list_x(il,ib),bc_list_y(il,ib))*
     .    hx(leftix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .       leftiy(bc_list_x(il,ib),bc_list_y(il,ib)))+
     .    fun(leftix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        leftiy(bc_list_x(il,ib),bc_list_y(il,ib)))*
     .    hx(bc_list_x(il,ib),bc_list_y(il,ib)))
     .  /(hx(bc_list_x(il,ib),bc_list_y(il,ib))+
     .    hx(leftix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .       leftiy(bc_list_x(il,ib),bc_list_y(il,ib))))
      elseif (bcchar(ib).eq.'W') then
        bcint =
     .   (fun(bc_list_x(il,ib),bc_list_y(il,ib))*
     .    hx(rightix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .       rightiy(bc_list_x(il,ib),bc_list_y(il,ib)))+
     .    fun(rightix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        rightiy(bc_list_x(il,ib),bc_list_y(il,ib)))*
     .    hx(bc_list_x(il,ib),bc_list_y(il,ib)))
     .  /(hx(bc_list_x(il,ib),bc_list_y(il,ib))+
     .    hx(rightix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .       rightiy(bc_list_x(il,ib),bc_list_y(il,ib))))
      else
        call xerrab ('bcint : Unknown orientation!')
      endif

* ..return
      return

*-----------------------------------------------------------------------
*.end bcint

      end function bcint
