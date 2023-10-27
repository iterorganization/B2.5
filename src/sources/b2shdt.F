*-----------------------------------------------------------------------
*.specification

      subroutine b2shdt (nCv, ns, dtim,                       !srv 07.07.21
     &  cvVol, na0, tep, ti0, tn0, ne0, ni0, nn0, kinrgy0,
     &  switch, mpg, shedt, shidt, shndt)
      use b2mod_types
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_balance !djm Jan2017
     1 , only : b2srdt_she0to3, b2srdt_shi0to3, balance_netcdf
      use b2mod_switches
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      type (switches), intent (in) :: switch
      type (mapping), intent (in) :: mpg
      integer nCv, ns
      real (kind=R8) ::
     *  dtim, cvVol(nCv), na0(nCv,0:ns-1),
     &  tep(nCv), ti0(nCv), tn0(nCv),
     &  ne0(nCv), ni0(nCv,0:1), nn0(nCv),
     &  kinrgy0(nCv,0:ns-1)
*   ..input/output arguments
      real (kind=R8) ::
     &  shedt(nCv,0:3), shidt(nCv,0:3), shndt(nCv,0:3)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SHDT modifies the source coefficients to account for the
*     implicit time discretisation for Te, Ti and Tn equations.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv, is, k
      character chk*1                                                !srv 11.09.09 }
      real (kind=R8) ::
     *   t0, ttim
*   ..procedures
      external xertst
      external b2xvsg
*   ..initialisation
      data ncall/0/                                                  !xpb

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2shdt')

*   ..test sign of dtim
      call xertst (0.0_R8.lt.dtim, 'faulty argument dtim')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test old state
       call b2xvsg (nCv*ns, na0, 1, 'na0', '.gt.')
       call b2xvsg (nCv, tep, 1, 'tep', '.gt.')
       call b2xvsg (nCv, ti0, 1, 'ti0', '.gt.')
       call b2xvsg (nCv, ne0, 1, 'ne0', '.gt.')
       call b2xvsg (2*nCv, ni0, 1, 'ni0', '.gt.')
       call b2xvsg (nCv*ns, kinrgy0, 1, 'kinrgy0', '.ge.')
      endif

*   ..modify shedt, shidt, shndt

      do iCv = 1, nCv
        ttim=dtim*dtee(mpg%cvReg(iCv))*time_factor(iCv)
        t0 = (switch%b2srdt_phm3/ttim)*cvVol(iCv)
        shedt(iCv,0) = 1.5_R8*t0*ne0(iCv)*tep(iCv)               !srv 06.01.21 {
        shedt(iCv,1) = 0.0_R8
        shedt(iCv,2) = 0.0_R8
        shedt(iCv,3) = -1.5_R8*t0                                    !srv 06.01.21 }
      enddo

! TODO energy equations !! tn_style !!
      do iCv = 1, nCv
        ttim=dtim*dtei(mpg%cvReg(iCv))*time_factor(iCv)
        t0 = (switch%b2srdt_phm3/ttim)*cvVol(iCv)
        shidt(iCv,0) = 1.5_R8*t0*ni0(iCv,0)*ti0(iCv)             !srv 06.01.21 {
        shidt(iCv,1) = 0.0_R8
        shidt(iCv,2) = 0.0_R8
        shidt(iCv,3) = -1.5_R8*t0                                    !srv 06.01.21 }
      enddo

      do iCv = 1, nCv
        ttim=dtim*dten(mpg%cvReg(iCv))*time_factor(iCv)
        t0 = (switch%b2srdt_phm3/ttim)*cvVol(iCv)
        shndt(iCv,0) = 1.5_R8*t0*nn0(iCv)*tn0(iCv)             !srv 06.01.21 {
        shndt(iCv,1) = 0.0_R8
        shndt(iCv,2) = 0.0_R8
        shndt(iCv,3) = -1.5_R8*t0                                 !srv 06.01.21 }
      enddo

  !! TODO_B2MNDT keep BORIS ??? otherwise remove kinrgy

cdjm Store the total change in linearised sources for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO        b2srdt_she0to3 = shedt
!WG_TODO        b2srdt_shi0to3 = shidt
!WG_TODO      endif

*    The diagnostic output of the time derivative will be done
*    in the subroutines where the solution is calculated.
      if (switch%b2shdt_iout.eq.1) then
       do k=0,3
        write (chk,'(i1)') k
        call my_out_us(70,nCv,shedt(:,k),'b2shdt_shedt'//chk)
        call my_out_us(70,nCv,shidt(:,k),'b2shdt_shidt'//chk)
        call my_out_us(70,nCv,shndt(:,k),'b2shdt_shndt'//chk)
       enddo
      endif                                                           !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2shdt

      end subroutine b2shdt
