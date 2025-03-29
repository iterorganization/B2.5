*-----------------------------------------------------------------------
*.specification

      subroutine b2srdt (nCv, ns, dtim, switch, geo, mpg, dv,
     &  na0, ua0, te0, ti0, tn0, ne0, ni0, nn0, kinrgy0, kt0, zt0,
     &  na, ua, te, ti, tn, ne, ni, nn, kinrgy, kt, zt,                 !srv 16.07.10
     &  sr, lout)              !srv 11.09.09
      use b2mod_types
      use b2mod_numerics_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_balance !djm Jan2017
     1 , only : b2srdt_sna0to1,b2srdt_smo0to3,b2srdt_she0to3,
     2          b2srdt_shi0to3, balance_netcdf
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2srdt
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      logical(1) :: lout                                                !srv 16.07.10
      integer nCv, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Derivatives), intent (in) :: dv
      real (kind=R8) ::
     *  dtim, na0(nCv,0:ns-1),
     &  ua0(nCv,0:ns-1), te0(nCv), ti0(nCv), tn0(nCv),
     &  ne0(nCv), ni0(nCv,0:1), nn0(nCv),
     &  kinrgy0(nCv,0:ns-1), kt0(nCv), zt0(nCv),
     &  na(nCv,0:ns-1), ua(nCv,0:ns-1),              !srv 16.07.10
     &  te(nCv), ti(nCv), tn(nCv), ne(nCv),
     &  ni(nCv,0:1), nn(nCv), kinrgy(nCv,0:ns-1), kt(nCv), zt(nCv)
*   ..input/output arguments
      type (B2Source), intent (inout) :: sr
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SRDT modifies the source coefficients to account for the
*     implicit time discretisation.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is, k
      character chns*3, chk*1                                        !srv 11.09.09 }
      real (kind=R8) :: t0, ttim
      real (kind=R8) :: wrk0(nCv), mult_dt(nCv)
      real (kind=R8) :: ttim2, src_tot, content, smt(0:3)
      real (kind=R8) :: ttim2_na(nCv,0:ns-1), ttim2_ua(nCv,0:ns-1),   ! only for iout = 1 + var_style = 2
     &                  ttim2_ti(nCv), ttim2_te(nCv),
     &                  ttim_na(nCv,0:ns-1), ttim_ua(nCv,0:ns-1),       
     &                  ttim_ti(nCv), ttim_te(nCv)
*   ..procedures
      external xertst
      external b2xvsg
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2srdt')
*   ..test nCv, ns
      call xertst (0.lt.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test sign of dtim
      call xertst (0.0_R8.lt.dtim, 'faulty argument dtim')
*   ..extensive tests on first few calls
      if (ncall_b2srdt.lt.3) then
*    ..test old state
       call b2xvsg (nCv*ns, na0, 1, 'na0', '.gt.')
       call b2xvsg (nCv, te0, 1, 'te0', '.gt.')
       call b2xvsg (nCv, ti0, 1, 'ti0', '.gt.')
       call b2xvsg (nCv, tn0, 1, 'tn0', '.gt.')
       call b2xvsg (nCv, ne0, 1, 'ne0', '.gt.')
       call b2xvsg (nCv, nn0, 1, 'nn0', '.gt.')
       call b2xvsg (2*nCv, ni0, 1, 'ni0', '.gt.')
       call b2xvsg (nCv*ns, kinrgy0, 1, 'kinrgy0', '.ge.')
       call b2xvsg (nCv, kt0, 1, 'kt0', '.ge.')
       call b2xvsg (nCv, zt0, 1, 'zt0', '.ge.')
       call xertst (maxval(abs(ua0)).lt.c, 'Supra-luminal velocity !')
*    ..test current state
       call b2xvsg (nCv*ns, na, 1, 'na', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, tn, 1, 'tn', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, nn, 1, 'nn', '.gt.')
       call b2xvsg (2*nCv, ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv*ns, kinrgy, 1, 'kinrgy', '.ge.')
       call b2xvsg (nCv, kt, 1, 'kt', '.ge.')
       call b2xvsg (nCv, zt, 1, 'zt', '.ge.')
       call xertst (maxval(abs(ua)).lt.c, 'Supra-luminal velocity !')
      endif

* ..Optionally, calculate a local (cell-based) multiplier to the timestep.
*   mult_dt will be 1.0 if b2mndt_var_style = 0 (default). 
      call get_mult_dt(switch%b2mndt_var_style, dtim, dv, geo, mpg, 
     &     nCv, ns,  te, switch%b2mndt_var_te_min,
     &     switch%b2mndt_var_te_max, switch%b2mndt_var_mult, 
     &     switch%b2mndt_var_style3_Cmax, mult_dt)

* ..include the contributions from dtim
*   ..modify snadt
      k = 1
      dtee_fts = 0.0_R8
      dtei_fts = 0.0_R8
      dtco_fts = 0.0_R8
      dtmo_fts = 0.0_R8
      do while (dtfts(k) .ne. 0 .and. k .le. DEF_NYD)
        dtee_fts(dtfts(k)) = dtee_ft(k)
        dtei_fts(dtfts(k)) = dtei_ft(k)
        dtco_fts(dtfts(k)) = dtco_ft(k)
        dtmo_fts(dtfts(k)) = dtmo_ft(k)
        k = k + 1
      enddo

#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$omp parallel
!$omp& default(none)
!$omp& private(is,iCv,ttim,t0, iFt, iFc)
!$omp& shared(nCv,ns,time_factor, geo, mpg, sr, lout)
!$omp& shared(dtim,dtco,dtmo,dtee,dtei, switch)
!$omp& shared(dtco_fts,dtmo_fts,dtee_fts,dtei_fts)
!$omp& shared(kinrgy0)
!$omp& shared(na0,ne0,ni0,ua0,ti0,te0,am, kt0, zt0, nn0, tn0)
!$omp& shared(ttim_na,ttim2_na,ttim_ua,ttim2_ua)
!$omp& shared(ttim_te,ttim2_te,ttim_ti,ttim2_ti,mult_dt)
!$omp& private(ttim2,content,src_tot)
#endif
      do is = 0, ns-1
        do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO           if(iftimbound(ix,iy) .ne. 0) cycle
!WG_TODO           if((guard_flows.eq.2 .or.
!WG_TODO     &        (is_neutral(is) .and. use_eirene.ne.0)) .and.
!WG_TODO     &       (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
!WG_TODO     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)) cycle
          ttim=dtim*dtco(is,mpg%cvReg(iCv))*time_factor(iCv)*
     &         mult_dt(iCv)
          if (switch%b2mndt_var_style.eq.2) then !! Timestep directly based on ratio of source to particle content
            src_tot = sr%sna(iCv,0,is)+sr%sna(iCv,1,is)*na0(iCv,is)
            content = geo%cvVol(iCv)*na0(iCv,is)
            ttim2 = content/abs(src_tot)/
     &              switch%b2mndt_var_style2_SF
            ttim  = min(max(ttim2,ttim),ttim*switch%b2mndt_var_mult)
            if (switch%b2srdt_iout.eq.1) then
              ttim2_na(iCv,is) = ttim2
              ttim_na(iCv,is) = ttim
            endif
          endif
          sr%snadt(iCv,0,is) =
     &     (switch%b2srdt_phm0/ttim)*geo%cvVol(iCv)*na0(iCv,is)
          sr%snadt(iCv,1,is) = -
     &     (switch%b2srdt_phm0/ttim)*geo%cvVol(iCv)
          sr%sna(iCv,:,is) = sr%sna(iCv,:,is)+sr%snadt(iCv,:,is)
        enddo
      enddo

*   ..modify smodt
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$omp DO
#endif
      do is = 0, ns-1
        do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO            if(iftimbound(ix,iy) .ne. 0) cycle

          ttim=dtim*dtmo(is,mpg%cvReg(iCv))*time_factor(iCv)*
     &         mult_dt(iCv)
          if (switch%b2mndt_var_style.eq.2) then !! Timestep directly based on ratio of source to parallel momentum content
            smt = sr%smq(iCv,0:3,is) + sr%smo(iCv,0:3,is)
            src_tot = smt(0) + 
     &                smt(1)*ua0(iCv,is) + 
     &                smt(2)*am(is)*mp*na0(iCv,is) +
     &                smt(3)*am(is)*mp*ua0(iCv,is)*na0(iCv,is)
            content = geo%cvVol(iCv)*geo%cvHz(iCv)*
     &                am(is)*mp*ua0(iCv,is)*na0(iCv,is)
            ttim2 = abs(content/src_tot)/
     &              switch%b2mndt_var_style2_SF
            ttim  = min(max(ttim2,ttim),ttim*switch%b2mndt_var_mult)
            if (switch%b2srdt_iout.eq.1) then
              ttim2_ua(iCv,is) = ttim2
              ttim_ua(iCv,is) = ttim
            endif
          endif
          t0 = (switch%b2srdt_phm1/ttim)*geo%cvVol(iCv)*geo%cvHz(iCv)   !srv 09.01.01
          sr%smodt(iCv,0,is) = t0*am(is)*mp*ua0(iCv,is)*na0(iCv,is)
          sr%smodt(iCv,1,is) = 0.0_R8
          sr%smodt(iCv,2,is) = 0.0_R8
          sr%smodt(iCv,3,is) = -t0
          sr%smo(iCv,:,is) = sr%smo(iCv,:,is)+sr%smodt(iCv,:,is)
        enddo
      enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$omp end do

*   ..modify shedt, shidt, shndt, sktdt, sztdt, schdt, snedt
#ifndef NO_OPENMP_B2SRDT
!$OMP DO
#endif
      do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO          if(iftimbound(ix,iy) .ne. 0) cycle

        ttim=dtim*dtee(mpg%cvReg(iCv))*time_factor(iCv)*
     &         mult_dt(iCv)
        if (switch%b2mndt_var_style.eq.2) then !! Timestep directly based on ratio of source to electron heat content
            src_tot = sr%she(iCv,0) + 
     &                sr%she(iCv,1)*te0(iCv) + 
     &                sr%she(iCv,2)*ne0(iCv) +
     &                sr%she(iCv,3)*te0(iCv)*ne0(iCv)
            content = 1.5_R8*ne0(iCv)*te0(iCv)*geo%cvVol(iCv)
            ttim2 = content/abs(src_tot)/
     &              switch%b2mndt_var_style2_SF
            ttim  = min(max(ttim2,ttim),ttim*switch%b2mndt_var_mult)
            if (switch%b2srdt_iout.eq.1) then
              ttim2_te(iCv) = ttim2
              ttim_te(iCv) = ttim
            endif
        endif
        t0 = (switch%b2srdt_phm3/ttim)*geo%cvVol(iCv)
        sr%shedt(iCv,0) = 1.5_R8*t0*ne0(iCv)*te0(iCv)
        sr%shedt(iCv,1) = 0.0_R8
        sr%shedt(iCv,2) = 0.0_R8
        sr%shedt(iCv,3) = -1.5_R8*t0
        sr%she(iCv,:) = sr%she(iCv,:)+sr%shedt(iCv,:)

        ttim=dtim*dtei(mpg%cvReg(iCv))*time_factor(iCv)*
     &         mult_dt(iCv)
        if (switch%b2mndt_var_style.eq.2) then !! For now keep timestep for kt and zt same is for ti
            src_tot = sr%shi(iCv,0) + 
     &                sr%shi(iCv,1)*ti0(iCv) + 
     &                sr%shi(iCv,2)*ni0(iCv,0) +
     &                sr%shi(iCv,3)*ti0(iCv)*ni0(iCv,0)
            content = 1.5_R8*ne0(iCv)*te0(iCv)*geo%cvVol(iCv)
            ttim2 = content/abs(src_tot)/
     &              switch%b2mndt_var_style2_SF
            ttim  = min(max(ttim2,ttim),ttim*switch%b2mndt_var_mult)
            if (switch%b2srdt_iout.eq.1) then
              ttim2_ti(iCv) = ttim2
              ttim_ti(iCv) = ttim
            endif
        endif
        t0 = (switch%b2srdt_phm3/ttim)*geo%cvVol(iCv)
        sr%sktdt(iCv,0) = t0*ni0(iCv,1)*kt0(iCv)
        sr%sktdt(iCv,1) = 0.0_R8
        sr%sktdt(iCv,2) = 0.0_R8
        sr%sktdt(iCv,3) = -t0
        sr%sktdt(iCv,:) = sr%skt(iCv,:)+sr%sktdt(iCv,:)

        sr%sztdt(iCv,0) = t0*ni0(iCv,1)*zt0(iCv)
        sr%sztdt(iCv,1) = 0.0_R8
        sr%sztdt(iCv,2) = 0.0_R8
        sr%sztdt(iCv,3) = -t0
        sr%sztdt(iCv,:) = sr%szt(iCv,:)+sr%sztdt(iCv,:)
      enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$omp end do

!$omp do
#endif
      do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO          if(iftimbound(ix,iy) .ne. 0) cycle

        ttim=dtim*dtei(mpg%cvReg(iCv))*time_factor(iCv)*
     &         mult_dt(iCv)
        if (switch%b2mndt_var_style.eq.2) then !! Timestep directly based on ratio of source to ion heat content
            src_tot = sr%shi(iCv,0) + 
     &                sr%shi(iCv,1)*ti0(iCv) + 
     &                sr%shi(iCv,2)*ni0(iCv,0) +
     &                sr%shi(iCv,3)*ti0(iCv)*ni0(iCv,0)
            content = 1.5_R8*ne0(iCv)*te0(iCv)*geo%cvVol(iCv)
            ttim2 = content/abs(src_tot)/
     &              switch%b2mndt_var_style2_SF
            ttim  = min(max(ttim2,ttim),ttim*switch%b2mndt_var_mult)
        endif
        t0 = (switch%b2srdt_phm3/ttim)*geo%cvVol(iCv)

        if (switch%tn_style.eq.0) then

          sr%shidt(iCv,0) = 1.5_R8*t0*ni0(iCv,0)*ti0(iCv)
          sr%shidt(iCv,1) = 0.0_R8
          sr%shidt(iCv,2) = 0.0_R8
          sr%shidt(iCv,3) = -1.5_R8*t0
          if (switch%BoRiS.eq.1.0_R8) then
            do is = 0, ns-1
              sr%shidt(iCv,0) = sr%shidt(iCv,0)+t0*
     &         kinrgy0(iCv,is)*na0(iCv,is)
              sr%shidt(iCv,2) = sr%shidt(iCv,2)-t0*
     &         kinrgy0(iCv,is)*na0(iCv,is)/
     &         ni0(iCv,0)
            enddo
          endif

        elseif (switch%tn_style.eq.1) then

          sr%shidt(iCv,0) = 1.5_R8*t0*ni0(iCv,0)*ti0(iCv)
          sr%shidt(iCv,1) = 0.0_R8
          sr%shidt(iCv,2) = 0.0_R8
          sr%shidt(iCv,3) = -1.5_R8*t0

        else

          sr%shidt(iCv,0) = 1.5_R8*t0*ni0(iCv,0)*ti0(iCv)
          sr%shidt(iCv,1) = 0.0_R8
          sr%shidt(iCv,2) = 0.0_R8
          sr%shidt(iCv,3) = -1.5_R8*t0

          sr%shndt(iCv,0) = 1.5_R8*t0*nn0(iCv)*tn0(iCv)
          sr%shndt(iCv,1) = 0.0_R8
          sr%shndt(iCv,2) = 0.0_R8
          sr%shndt(iCv,3) = -1.5_R8*t0
          sr%shn(iCv,:) = sr%shn(iCv,:)+sr%shndt(iCv,:)
        endif

        sr%shi(iCv,:) = sr%shi(iCv,:)+sr%shidt(iCv,:)
      enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$omp end do

!$omp do
#endif
      do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO          if(iftimbound(ix,iy) .ne. 0) cycle


cxpb  SCH should normally remain unchanged (phm4.eq.0.0)
cxpb  This treatment added as a fix to get the ExB runs going
        t0 = (switch%b2srdt_phm4/dtim)
     &               *geo%cvVol(iCv)*qe*time_factor(iCv)                !xpb
        sr%schdt(iCv,0) = t0*(ni0(iCv,0)-ne0(iCv))                      !xpb
        sr%schdt(iCv,1) = 0.0_R8
        sr%schdt(iCv,2) = 0.0_R8
        sr%schdt(iCv,3) = -t0                                           !xpb
        sr%sch(iCv,:) = sr%sch(iCv,:)+sr%schdt(iCv,:)
      enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$omp end do

!$omp do
#endif
      do iCv = 1, nCv
!WG_TODO          if(iftimbound(ix,iy) .ne. 0) cycle
        t0 = (switch%b2srdt_phm5/dtim)*geo%cvVol(iCv)
     &                                  *time_factor(iCv)               !xpb
        sr%snedt(iCv,0) = t0*ne0(iCv)
        sr%snedt(iCv,1) = -t0
        sr%sne(iCv,:) = sr%sne(iCv,:)+sr%snedt(iCv,:)
      enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$omp end do
!$omp end parallel
#endif

cdjm Store the total change in linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2srdt_sna0to1 = sr%snadt
        b2srdt_smo0to3 = sr%smodt
        b2srdt_she0to3 = sr%shedt
        b2srdt_shi0to3 = sr%shidt
      endif

      if ((switch%b2npmo_iout.eq.1.or.
     &                   switch%iout_b2wdat .eq. 4).and.lout) then      !srv 11.09.09 {
       do is = 0, ns-1
         write (chns,'(i3.3)') is
         wrk0=sr%smodt(:,0,is)+
     &        sr%smodt(:,1,is)*ua(:,is)+
     &        sr%smodt(:,2,is)*na(:,is)*mp*am(is)+
     &        sr%smodt(:,3,is)*ua(:,is)*na(:,is)*mp*am(is)
         call my_out_us(70,nCv,0,wrk0,'b2srdt_smodt'//chns)
       enddo
      endif
      if ((switch%b2npco_iout.eq.1
     &    .or. switch%iout_b2wdat.eq.4) .and. lout) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          wrk0=sr%snadt(:,0,is)+sr%snadt(:,1,is)*na(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2srdt_snadt'//chns)
        enddo
      endif
      if ((switch%b2npht_iout.eq.1 .or.
     &                      switch%iout_b2wdat.eq.4) .and. lout) then
        wrk0=sr%shedt(:,0)+sr%shedt(:,1)*te+
     &       sr%shedt(:,2)*ne + sr%shedt(:,3)*te*ne
        call my_out_us(70,nCv,0,wrk0,'b2srdt_shedt')
        wrk0=sr%shidt(:,0)+sr%shidt(:,1)*ti+
     &       sr%shidt(:,2)*ni(:,0)+
     &       sr%shidt(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,0,wrk0,'b2srdt_shidt')
      endif

      if (switch%b2srdt_iout.eq.1 .and. lout) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        do k=0,3
         write (chk,'(i1)') k
         call my_out_us(70,nCv,0,sr%smodt(1,k,is),
     &                       'b2srdt_smodt'//chk//chns)
        enddo
        do k=0,1
         write (chk,'(i1)') k
         call my_out_us(70,nCv,0,sr%snadt(1,k,is),
     &                       'b2srdt_snadt'//chk//chns)
        enddo
       enddo
       do k=0,3
        write (chk,'(i1)') k
        call my_out_us(70,nCv,0,sr%shedt(1,k),'b2srdt_shedt'//chk)
        call my_out_us(70,nCv,0,sr%shidt(1,k),'b2srdt_shidt'//chk)
        if (switch%tn_style.ge.2)
     &  call my_out_us(70,nCv,0,sr%shndt(1,k),'b2srdt_shndt'//chk)
        call my_out_us(70,nCv,0,sr%schdt(1,k),'b2srdt_schdt'//chk)
        call my_out_us(70,nCv,0,sr%sktdt(1,k),'b2srdt_sktdt'//chk)
        call my_out_us(70,nCv,0,sr%sztdt(1,k),'b2srdt_sztdt'//chk)
       enddo
      endif                                                           !srv 11.09.09 }

      if (switch%b2srdt_iout.eq.1.and.switch%b2mndt_var_style.eq.2) then
c print out the variable local timesteps for b2mndt_var_style = 2.
c ttim2 = unlimited value
c ttim = value limited by minumum timestep (=baseline timestep) and by maximum timestep (=baseline timestep x multiplier)
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          wrk0 = ttim2_na(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2srdt_ttim2_na'//chns)
          wrk0 = ttim_na(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2srdt_ttim_na'//chns)
          wrk0 = ttim2_ua(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2srdt_ttim2_ua'//chns)
          wrk0 = ttim_ua(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2srdt_ttim_ua'//chns)
        enddo
        call my_out_us(70,nCv,0,ttim2_te,'b2srdt_ttim2_te')
        call my_out_us(70,nCv,0,ttim_te,'b2srdt_ttim_te')
        call my_out_us(70,nCv,0,ttim2_ti,'b2srdt_ttim2_ti')
        call my_out_us(70,nCv,0,ttim_ti,'b2srdt_ttim_ti')
      endif
      if (switch%b2srdt_iout.eq.1.and.switch%b2mndt_var_style.eq.1) then
        call my_out_us(70,nCv,0,mult_dt,'b2srdt_mult_dt')
      endif

* ..return
      ncall_b2srdt = ncall_b2srdt+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srdt

      end subroutine b2srdt
