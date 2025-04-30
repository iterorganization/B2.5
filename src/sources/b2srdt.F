*-----------------------------------------------------------------------
*.specification

      subroutine b2srdt (nCv, ns, dtim, switch, geo, mpg,
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
#ifndef FIXED_POINT
      use b2mod_balance !djm Jan2017
     1 , only : b2srdt_sna0to1,b2srdt_smo0to3,b2srdt_she0to3,
     2          b2srdt_shi0to3, balance_netcdf
#endif
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
      integer iCv, is, k, iFt, iFc
      character chns*3, chk*1                                        !srv 11.09.09 }
      real (kind=R8) :: t0, ttim
      real (kind=R8) :: wrk0(nCv)
      real (kind=R8) :: dtee_fts(mpg%nFt), dtei_fts(mpg%nFt),
     &                  dtco_fts(mpg%nFt), dtmo_fts(mpg%nFt)
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
       call xertst (maxval(dabs(ua0)).lt.c, 'Supra-luminal velocity !')
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
       call xertst (maxval(dabs(ua)).lt.c, 'Supra-luminal velocity !')
      endif

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
      do is = 0, ns-1
        do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO           if(iftimbound(ix,iy) .ne. 0) cycle
!WG_TODO           if((guard_flows.eq.2 .or.
!WG_TODO     &        (is_neutral(is) .and. use_eirene.ne.0)) .and.
!WG_TODO     &       (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
!WG_TODO     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)) cycle
          ttim=dtim*dtco(is,mpg%cvReg(iCv))*time_factor(iCv)
          iFt = mpg%cvFt(iCv)
          if (iFt.eq.0 .and. iCv.gt.mpg%nCi) then
            iFc = mpg%cvFc(mpg%cvFcP(iCv,1))
            iFt = mpg%cvFt(minval(mpg%fcCv(iFc,:)))
          endif
          if (iFt.ne.0) then
            if (dtco_fts(iFt).gt.0.0_R8) then
              ttim=dtim*dtco_fts(iFt)*time_factor(iCv)
            endif
          endif
          sr%snadt(iCv,0,is) =
     &     (switch%b2srdt_phm0/ttim)*geo%cvVol(iCv)*na0(iCv,is)
          sr%snadt(iCv,1,is) = -
     &     (switch%b2srdt_phm0/ttim)*geo%cvVol(iCv)
          if (.not.lout)
     &     sr%sna(iCv,:,is) = sr%sna(iCv,:,is)+sr%snadt(iCv,:,is)
        enddo
      enddo

*   ..modify smodt
      do is = 0, ns-1
        do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO            if(iftimbound(ix,iy) .ne. 0) cycle

          ttim=dtim*dtmo(is,mpg%cvReg(iCv))*time_factor(iCv)
          iFt = mpg%cvFt(iCv)
          if (iFt.eq.0 .and. iCv.gt.mpg%nCi) then
            iFc = mpg%cvFc(mpg%cvFcP(iCv,1))
            iFt = mpg%cvFt(minval(mpg%fcCv(iFc,:)))
          endif
          if (iFt.ne.0) then
            if (dtmo_fts(iFt).gt.0.0_R8) then
              ttim=dtim*dtmo_fts(iFt)*time_factor(iCv)
            endif
          endif
          t0 = (switch%b2srdt_phm1/ttim)*geo%cvVol(iCv)*geo%cvHz(iCv)   !srv 09.01.01
          sr%smodt(iCv,0,is) = t0*am(is)*mp*ua0(iCv,is)*na0(iCv,is)
          sr%smodt(iCv,1,is) = 0.0_R8
          sr%smodt(iCv,2,is) = 0.0_R8
          sr%smodt(iCv,3,is) = -t0
          if (.not.lout)
     &     sr%smo(iCv,:,is) = sr%smo(iCv,:,is)+sr%smodt(iCv,:,is)
        enddo
      enddo

*   ..modify shedt, shidt, shndt, sktdt, sztdt, schdt, snedt
      do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO          if(iftimbound(ix,iy) .ne. 0) cycle

        ttim=dtim*dtee(mpg%cvReg(iCv))*time_factor(iCv)
        iFt = mpg%cvFt(iCv)
        if (iFt.eq.0 .and. iCv.gt.mpg%nCi) then
          iFc = mpg%cvFc(mpg%cvFcP(iCv,1))
          iFt = mpg%cvFt(minval(mpg%fcCv(iFc,:)))
       endif
        if (iFt.ne.0) then
          if (dtee_fts(iFt).gt.0.0_R8) then
            ttim=dtim*dtee_fts(iFt)*time_factor(iCv)
          endif
        endif
        t0 = (switch%b2srdt_phm3/ttim)*geo%cvVol(iCv)
        sr%shedt(iCv,0) = 1.5_R8*t0*ne0(iCv)*te0(iCv)
        sr%shedt(iCv,1) = 0.0_R8
        sr%shedt(iCv,2) = 0.0_R8
        sr%shedt(iCv,3) = -1.5_R8*t0
        if (.not.lout)
     &   sr%she(iCv,:) = sr%she(iCv,:)+sr%shedt(iCv,:)

        ttim=dtim*dtei(mpg%cvReg(iCv))*time_factor(iCv)
        if (iFt.ne.0) then
          if (dtei_fts(iFt) .gt. 0.0_R8) then
            ttim=dtim*dtei_fts(iFt)*time_factor(iCv)
          endif
        endif
        t0 = (switch%b2srdt_phm3/ttim)*geo%cvVol(iCv)
        sr%sktdt(iCv,0) = t0*ni0(iCv,1)*kt0(iCv)
        sr%sktdt(iCv,1) = 0.0_R8
        sr%sktdt(iCv,2) = 0.0_R8
        sr%sktdt(iCv,3) = -t0
        if (.not.lout)
     &   sr%skt(iCv,:) = sr%skt(iCv,:)+sr%sktdt(iCv,:)

        sr%sztdt(iCv,0) = t0*ni0(iCv,1)*zt0(iCv)
        sr%sztdt(iCv,1) = 0.0_R8
        sr%sztdt(iCv,2) = 0.0_R8
        sr%sztdt(iCv,3) = -t0
        if (.not.lout)
     &   sr%szt(iCv,:) = sr%szt(iCv,:)+sr%sztdt(iCv,:)
      enddo

      do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO          if(iftimbound(ix,iy) .ne. 0) cycle

        ttim=dtim*dtei(mpg%cvReg(iCv))*time_factor(iCv)
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
          if (.not.lout)
     &     sr%shn(iCv,:) = sr%shn(iCv,:)+sr%shndt(iCv,:)
        endif

        if (.not.lout)
     &   sr%shi(iCv,:) = sr%shi(iCv,:)+sr%shidt(iCv,:)
      enddo

      do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO          if(iftimbound(iCv) .ne. 0) cycle
cxpb  SCH should normally remain unchanged (phm4.eq.0.0)
cxpb  This treatment added as a fix to get the ExB runs going
        t0 = (switch%b2srdt_phm4/dtim)
     &               *geo%cvVol(iCv)*qe*time_factor(iCv)                !xpb
        sr%schdt(iCv,0) = t0*(ni0(iCv,0)-ne0(iCv))                      !xpb
        sr%schdt(iCv,1) = 0.0_R8
        sr%schdt(iCv,2) = 0.0_R8
        sr%schdt(iCv,3) = -t0                                           !xpb
        if (.not.lout)
     &   sr%sch(iCv,:) = sr%sch(iCv,:)+sr%schdt(iCv,:)
      enddo

      do iCv = 1, nCv
!WG_TODO          if(iftimbound(iCv) .ne. 0) cycle
        t0 = (switch%b2srdt_phm5/dtim)*geo%cvVol(iCv)
     &                                  *time_factor(iCv)               !xpb
        sr%snedt(iCv,0) = t0*ne0(iCv)
        sr%snedt(iCv,1) = -t0
        if (.not.lout)
     &   sr%sne(iCv,:) = sr%sne(iCv,:)+sr%snedt(iCv,:)
      enddo

#ifndef FIXED_POINT
cdjm Store the total change in linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2srdt_sna0to1 = sr%snadt
        b2srdt_smo0to3 = sr%smodt
        b2srdt_she0to3 = sr%shedt
        b2srdt_shi0to3 = sr%shidt
      endif
#endif

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

* ..return
      ncall_b2srdt = ncall_b2srdt+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srdt

      end subroutine b2srdt
