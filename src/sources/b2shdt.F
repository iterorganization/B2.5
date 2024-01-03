*-----------------------------------------------------------------------
*.specification

      subroutine b2shdt (nCv, nCi, ns, dtim, BoRiS,                    !srv 07.07.21
     &  cvVol, nap, tep, tip, tnp, ktp, ztp, ne, ni, nn, kinrgyp,
     &  switch, mpg, shedt, shidt, shndt, sktdt, sztdt)
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
      integer nCv, nCi, ns, BoRiS
      real (kind=R8) ::
     *  dtim, cvVol(nCv), nap(nCv,0:ns-1),
     &  tep(nCv), tip(nCv), tnp(nCv),
     &  ktp(nCv), ztp(nCv),
     &  ne(nCv), ni(nCv,0:1), nn(nCv),
     &  kinrgyp(nCv,0:ns-1)
*   ..input/output arguments
      real (kind=R8) ::
     &  shedt(nCv,0:3), shidt(nCv,0:3), shndt(nCv,0:3),
     &  sktdt(nCv,0:3), sztdt(nCv,0:3)
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
       call b2xvsg (nCv*ns, nap, 1, 'nap', '.gt.')
       call b2xvsg (nCv, tep, 1, 'tep', '.gt.')
       call b2xvsg (nCv, tip, 1, 'tip', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (2*nCv, ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, nn, 1, 'nn', '.gt.')
       call b2xvsg (nCv*ns, kinrgyp, 1, 'kinrgyp', '.ge.')
      endif
      if (switch%b2srdt_phm3.ne.1.0_R8) then 
        write(*,*)
     & 'Warning: b2srdt_phm3 should be 1.0 for time-dependent mode'
      endif

*   ..modify shedt, shidt, shndt

      do iCv = 1, nCi ! only internal cells
        ttim=dtim*dtee(mpg%cvReg(iCv))*time_factor(iCv)
        t0 = (switch%b2srdt_phm3/ttim)*cvVol(iCv)
        shedt(iCv,0) = 1.5_R8*t0*ne(iCv)*tep(iCv)               !srv 06.01.21 {
        shedt(iCv,1) = 0.0_R8
        shedt(iCv,2) = 0.0_R8
        shedt(iCv,3) = -1.5_R8*t0                               !srv 06.01.21 }
      enddo

      do iCv = 1, nCi ! only internal cells
        ttim=dtim*dtei(mpg%cvReg(iCv))*time_factor(iCv)
        t0 = (switch%b2srdt_phm3/ttim)*cvVol(iCv)
        if (switch%tn_style.eq.0) then
          shidt(iCv,0) = 1.5_R8*t0*ni(iCv,0)*tip(iCv)           !srv 06.01.21 {
          shidt(iCv,1) = 0.0_R8
          shidt(iCv,2) = 0.0_R8
          shidt(iCv,3) = -1.5_R8*t0                             !srv 06.01.21 }
          if (BoRiS.eq.1.0_R8) then
            do is = 0, ns-1
              shidt(iCv,0) = shidt(iCv,0)+t0*
     &           (kinrgyp(iCv,is)+tip(iCv))*nap(iCv,is)         !srv 10.10.22 {
              shidt(iCv,2) = shidt(iCv,2)-t0*
     &           (kinrgyp(iCv,is)+tip(iCv))*nap(iCv,is)/
     &           ni(iCv,0)                                      !srv 10.10.22 }
            enddo
          endif
        else
          shidt(iCv,0) = 1.5_R8*t0*ni(iCv,1)*tip(iCv)
          shidt(iCv,1) = 0.0_R8
          shidt(iCv,2) = 0.0_R8
          shidt(iCv,3) = -1.5_R8*t0*ni(iCv,1)/ni(iCv,0)          ! this is needed because the source scaling for shi remains based on ni(:,0), even though physically it involves ni(:,1)
                                                                 ! when using a pure ion or separate ion-neutral energy equation
          if (BoRiS.eq.1.0_R8) then
            call xertst('BoRiS not allowed with tn_style.ne.0')
          endif
        endif
      enddo

      if (switch%tn_style.eq.2) then
        do iCv = 1, nCi ! only internal cells
          ttim=dtim*dten(mpg%cvReg(iCv))*time_factor(iCv)
          t0 = (switch%b2srdt_phm3/ttim)*cvVol(iCv)
          shndt(iCv,0) = 1.5_R8*t0*nn(iCv)*tnp(iCv)               !srv 06.01.21 {
          shndt(iCv,1) = 0.0_R8
          shndt(iCv,2) = 0.0_R8
          shndt(iCv,3) = -1.5_R8*t0                               !srv 06.01.21 }
        enddo
      endif

      if (switch%solve_keps.gt.0) then
        do iCv = 1, nCi ! only internal cells
          ttim=dtim*time_factor(iCv)
          t0 = (switch%b2srdt_phm3/ttim)*cvVol(iCv)
          sktdt(iCv,0) = 1.5_R8*t0*ni(iCv,1)*ktp(iCv)
          sktdt(iCv,1) = 0.0_R8
          sktdt(iCv,2) = 0.0_R8
          sktdt(iCv,3) = -1.5_R8*t0
        enddo
      endif

      if (switch%solve_keps.gt.1) then
        do iCv = 1, nCi ! only internal cells
          ttim=dtim*time_factor(iCv)
          t0 = (switch%b2srdt_phm3/ttim)*cvVol(iCv)
          sztdt(iCv,0) = 1.5_R8*t0*ni(iCv,1)*ztp(iCv)
          sztdt(iCv,1) = 0.0_R8
          sztdt(iCv,2) = 0.0_R8
          sztdt(iCv,3) = -1.5_R8*t0
        enddo
      endif


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
        call my_out_us(70,nCv,sktdt(:,k),'b2shdt_sktdt'//chk)
        call my_out_us(70,nCv,sztdt(:,k),'b2shdt_sztdt'//chk)
       enddo
      endif                                                           !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2shdt

      end subroutine b2shdt
