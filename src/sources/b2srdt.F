*-----------------------------------------------------------------------
*.specification

      subroutine b2srdt (nCv, ns, dtim, BoRiS, geo, mpg,
     &  na0, ua0, te0, ti0, tn0, ne0, ni0, kinrgy0, kt0,
     &  na, ua, te, ti, tn, ne, ni, kinrgy, kt,                           !srv 16.07.10
     &  sr, lout)              !srv 11.09.09
      use b2mod_types 
      use b2mod_numerics_namelist 
      use b2mod_constants 
      use b2mod_b2cmpa 
      use b2us_geo
      use b2us_map
      use b2us_plasma
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     1 , only : b2srdt_sna0to1,b2srdt_smo0to3,b2srdt_she0to3,
!WG_TODO     2          b2srdt_shi0to3, balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      logical(1) :: lout                                             !srv 16.07.10
      integer nCv, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  dtim, BoRiS, na0(nCv,0:ns-1),
     &  ua0(nCv,0:ns-1), te0(nCv), ti0(nCv), tn0(nCv),
     &  ne0(nCv), ni0(nCv,0:1),
     &  kinrgy0(nCv,0:ns-1), kt0(nCv),
     &  na(nCv,0:ns-1), ua(nCv,0:ns-1),              !srv 16.07.10
     &  te(nCv), ti(nCv), tn(nCv), ne(nCv),
     &  ni(nCv,0:1), kinrgy(nCv,0:ns-1), kt(nCv)
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
      integer ncall, iCv, is, k, use_eirene
      integer, save :: iout = 0, iout_b2npc9 = 0, iout_b2npco = 0,   !srv 11.09.09 {
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npmo = 0, iout_b2npp7 = 0,
     &                 iout_b2nppo = 0, iout_b2wdat = 0,
     &                 tn_style = 0
      character chns*3, chk*1                                        !srv 11.09.09 }
      real (kind=R8) ::
     *   phm0, phm1, phm3, phm4, phm5, t0, ttim
      real (kind=R8) :: wrk0(nCv)
      save ncall, phm0, phm1, phm3, phm4, phm5, use_eirene           !xpb
      real (kind=R8) :: !djm Jan2017
     &  snatmp(nCv,0:1,0:ns-1), smotmp(nCv,0:3,0:ns-1), 
     &  shetmp(nCv,0:3), shitmp(nCv,0:3), shntmp(nCv,0:3)
*   ..procedures
      external subini, subend, xertst, ipgetr
      external b2xvsg
*   ..initialisation
      data ncall/0/, phm0/1.0_R8/, phm1/1.0_R8/, phm3/1.0_R8/,
     1  phm4/0.0_R8/, phm5/1.0_R8/, use_eirene/0/      		     !xpb

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2srdt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgetr ('b2srdt_phm0', phm0)
        call ipgetr ('b2srdt_phm1', phm1)
        call ipgetr ('b2srdt_phm3', phm3)
        call ipgetr ('b2srdt_phm4', phm4)			     !xpb
        call ipgetr ('b2srdt_phm5', phm5)		             !xpb
        call xertst (0.0_R8.le.phm0.and.0.0_R8.le.phm1		     !xpb
     &   .and.0.0_R8.le.phm3.and.0.0_R8.le.phm4.and.0.0_R8.le.phm5,  !xpb
     &	 'faulty internal parameter phm0, phm1, phm3, phm4, phm5')   !xpb
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2srdt_iout', iout)                            !srv 30.06.08
        call ipgeti ('b2npc9_iout', iout_b2npc9)                     !srv 11.09.09 {
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)                     !srv 11.09.09 }
        call ipgeti ('b2wdat_iout', iout_b2wdat)                     !iys 17.04.17
        call ipgeti ('b2mn_tn_style', tn_style)                  
      endif
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test sign of dtim
      call xertst (0.0_R8.lt.dtim, 'faulty argument dtim')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test old state
       call b2xvsg (nCv*ns, na0, 1, 'na0', '.gt.')
       call b2xvsg (nCv, te0, 1, 'te0', '.gt.')
       call b2xvsg (nCv, ti0, 1, 'ti0', '.gt.')
       call b2xvsg (nCv, ne0, 1, 'ne0', '.gt.')
       call b2xvsg (2*nCv, ni0, 1, 'ni0', '.gt.')
       call b2xvsg (nCv*ns, kinrgy0, 1, 'kinrgy0', '.ge.')
       call b2xvsg (nCv, kt0, 1, 'kt0', '.ge.')
      endif

cdjm Jan2017 Keep sources at the start to calculate total change for balance routines
      snatmp = sr%sna
      smotmp = sr%smo
      shetmp = sr%she
      shitmp = sr%shi
!!      shntmp = sr%shn

* ..include the contributions from dtim
*   ..modify snadt
      do is = 0, ns-1
        do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO           if(iftimbound(ix,iy) .ne. 0) cycle
!WG_TODO           if((guard_flows.eq.2 .or.
!WG_TODO     &        (is_neutral(is) .and. use_eirene.ne.0)) .and. 
!WG_TODO     &       (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
!WG_TODO     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)) cycle
         ttim=dtim*dtco(is,mpg%cvReg(iCv))*time_factor(iCv)
         sr%sna(iCv,0,is) = sr%sna(iCv,0,is)+
     &    (phm0/ttim)*geo%cvVol(iCv)*na0(iCv,is)
         sr%sna(iCv,1,is) = sr%sna(iCv,1,is)-
     &    (phm0/ttim)*geo%cvVol(iCv)
        enddo
      enddo
*   ..modify smodt
      do is = 0, ns-1
        do iCv = 1, nCv
*        smodt(-1,iy,0:3,is) = (unchanged at 0.0)
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO            if(iftimbound(ix,iy) .ne. 0) cycle

            ttim=dtim*dtmo(is,mpg%cvReg(iCv))*time_factor(iCv)
            t0 = (phm1/ttim)*geo%cvVol(iCv)*geo%cvHz(iCv)                    !srv 09.01.01
            sr%smo(iCv,0,is) = sr%smo(iCv,0,is)+
     1       t0*am(is)*mp*ua0(iCv,is)*
     &       na0(iCv,is)
*        smodt(iCv,1,is) = (unchanged)
*        smodt(iCv,2,is) = (unchanged)
            sr%smo(iCv,3,is) = sr%smo(iCv,3,is)-t0
        enddo
      enddo
*   ..modify shedt, shidt, schdt, snedt
      do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO          if(iftimbound(ix,iy) .ne. 0) cycle

          ttim=dtim*dtee(mpg%cvReg(iCv))*time_factor(iCv)
          t0 = (phm3/ttim)*geo%cvVol(iCv)
          sr%she(iCv,0) = sr%she(iCv,0)+
     &     1.5_R8*t0*ne0(iCv)*te0(iCv)
*         shedt(iCv,1) = (unchanged)
*         shedt(iCv,2) = (unchanged)
          sr%she(iCv,3) = sr%she(iCv,3)-1.5_R8*t0

          ttim=dtim*dtei(mpg%cvReg(iCv))*time_factor(iCv)
          t0 = (phm3/ttim)*geo%cvVol(iCv)
          sr%skt(iCv,0) = sr%skt(iCv,0)+
     &     t0*ni0(iCv,1)*kt0(iCv)
*         sktdt(iCv,1) = (unchanged)
*         sktdt(iCv,2) = (unchanged)
          sr%skt(iCv,3) = sr%skt(iCv,3)-t0
      enddo
      do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO          if(iftimbound(ix,iy) .ne. 0) cycle

          ttim=dtim*dtei(mpg%cvReg(iCv))*time_factor(iCv)
          t0 = (phm3/ttim)*geo%cvVol(iCv)
          if (tn_style.eq.0) then
            sr%shi(iCv,0) = sr%shi(iCv,0)+
     &       1.5_R8*t0*ni0(iCv,0)*ti0(iCv)
          elseif (tn_style.eq.1) then
            sr%shi(iCv,0) = sr%shi(iCv,0)+
     &       1.5_R8*t0*ni0(iCv,1)*ti0(iCv)
          else
            sr%shi(iCv,0) = sr%shi(iCv,0)+
     &       1.5_R8*t0*ni0(iCv,1)*ti0(iCv)
            sr%shn(iCv,0) = sr%shn(iCv,0)+
     &       1.5_R8*t0*(ni0(iCv,0)-ni0(iCv,1))*tn0(iCv)
          endif
          if (BoRiS.eq.1.0_R8) then
            if (tn_style.ne.0) then
              call xerrab('BoRiS only with default energy eq.')
            endif
            do is = 0, ns-1
              sr%shi(iCv,0) = sr%shi(iCv,0)+t0*
     &         kinrgy0(iCv,is)*na0(iCv,is)
              sr%shi(iCv,2) = sr%shi(iCv,2)-t0*
     &         kinrgy0(iCv,is)*na0(iCv,is)/
     &         ni0(iCv,0)
            enddo
          endif
*         shidt(iCv,1) = (unchanged)
*         shidt(iCv,2) = (otherwise unchanged)
          sr%shi(iCv,3) = sr%shi(iCv,3)-1.5_R8*t0

*         shndt(iCv,1) = (unchanged)
*         shndt(iCv,2) = (otherwise unchanged)
          sr%shn(iCv,3) = sr%shn(iCv,3)-1.5_R8*t0
      enddo
      do iCv = 1, nCv
csw 24feb2014 check for guard cell (SOLPS4)
!WG_TODO          if(iftimbound(ix,iy) .ne. 0) cycle


cxpb  SCH should normally remain unchanged (phm4.eq.0.0)
cxpb  This treatment added as a fix to get the ExB runs going
          t0 = (phm4/dtim)*geo%cvVol(iCv)*qe*time_factor(iCv)		 !xpb
          sr%sch(iCv,0) = sr%sch(iCv,0)+t0*(ni0(iCv,0)-ne0(iCv))	 !xpb
*         schdt(iCv,1) = (unchanged)
*         schdt(iCv,2) = (unchanged)
          sr%sch(iCv,3) = sr%sch(iCv,3)-t0				 !xpb
      enddo
      do iCv = 1, nCv
!WG_TODO          if(iftimbound(ix,iy) .ne. 0) cycle
          t0 = (phm5/dtim)*geo%cvVol(iCv)*time_factor(iCv)		 !xpb
          sr%sne(iCv,0) = sr%sne(iCv,0)+t0*ne0(iCv)
          sr%sne(iCv,1) = sr%sne(iCv,1)-t0
      enddo

cdjm Store the total change in linearised sources for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO        b2srdt_sna0to1 = snadt-snatmp
!WG_TODO        b2srdt_smo0to3 = smodt-smotmp
!WG_TODO        b2srdt_she0to3 = shedt-shetmp
!WG_TODO        b2srdt_shi0to3 = shidt-shitmp
!WG_TODO      endif

      if ((iout_b2npmo.eq.1 .or. iout_b2wdat .eq. 4) .and. lout) then                                     !srv 11.09.09 {
       do is = 0, ns-1 
         write (chns,'(i3.3)') is
         wrk0=sr%smo(:,0,is)+
     &        sr%smo(:,1,is)*ua(:,is)+
     &        sr%smo(:,2,is)*na(:,is)*mp*am(is)+
     &        sr%smo(:,3,is)*ua(:,is)*na(:,is)*mp*am(is)
         call my_out_us(70,nCv,0,wrk0,'b2srdt_smodt'//chns)
       enddo
      endif
      if ((iout_b2npc9.eq.1 .or. iout_b2npco.eq.1 .or. iout_b2wdat.eq.4)
     &             .and. lout) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          wrk0=sr%sna(:,0,is)+sr%sna(:,1,is)*na(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2srdt_snadt'//chns)
        enddo
      endif
      if ((iout_b2nph9.eq.1 .or. iout_b2npht.eq.1 .or. iout_b2wdat.eq.4)
     &             .and. lout) then
        wrk0=sr%she(:,0)+sr%she(:,1)*te+
     &       sr%she(:,2)*ne + sr%she(:,3)*te*ne
        call my_out_us(70,nCv,0,wrk0,'b2srdt_shedt')
        wrk0=sr%shi(:,0)+sr%shi(:,1)*ti+
     &       sr%shi(:,2)*ni(:,0)+
     &       sr%shi(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,0,wrk0,'b2srdt_shidt')
      endif

      if (iout.eq.1 .and. lout) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        do k=0,3
         write (chk,'(i1)') k
         call my_out_us(70,nCv,0,sr%smo(1,k,is),
     &                       'b2srdt_smodt'//chk//chns)
        enddo
        do k=0,1
         write (chk,'(i1)') k
         call my_out_us(70,nCv,0,sr%sna(1,k,is),
     &                       'b2srdt_snadt'//chk//chns)
        enddo
       enddo
       do k=0,3
        write (chk,'(i1)') k
        call my_out_us(70,nCv,0,sr%she(1,k),'b2srdt_shedt'//chk)
        call my_out_us(70,nCv,0,sr%shi(1,k),'b2srdt_shidt'//chk)
        call my_out_us(70,nCv,0,sr%sch(1,k),'b2srdt_schdt'//chk)
        call my_out_us(70,nCv,0,sr%skt(1,k),'b2srdt_sktdt'//chk)
       enddo
      endif                                                           !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srdt

      end subroutine b2srdt
