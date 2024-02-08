*-----------------------------------------------------------------------
*.specification

      subroutine b2scdt (nCv, nCi, isb, dtim,                        !srv 07.07.21
     &  cvVol, nap, switch, mpg, snadt)
      use b2mod_types
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_constants
      use b2mod_b2cmpa
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     1 , only : b2srdt_sna0to1, balance_netcdf
      use b2mod_switches
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      type (switches), intent (in) :: switch
      type (mapping), intent (in) :: mpg
      integer nCv, nCi, isb                                          !srv 01.07.08
      real (kind=R8) ::
     *  dtim, cvVol(nCv), nap(nCv)
*   ..input/output arguments
      real (kind=R8) ::
     *  snadt(nCv,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SCDT modifies the source coefficients to account for the
*     implicit time discretisation for the isb-th density equation.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv, k
      character chns*3, chk*1                                        !srv 11.09.09 }
      real (kind=R8) ::
     *   t0, ttim
      save ncall                                                     !xpb
*   ..procedures
      external xertst
      external b2xvsg
*   ..initialisation
      data ncall/0/                                                  !xpb

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2scdt')
*   ..set internal parameters on first call
      call xertst (0.le.isb, 'faulty argument isb')
*   ..test sign of dtim
      call xertst (0.0_R8.lt.dtim, 'faulty argument dtim')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test old state
       call b2xvsg (nCv, nap, 1, 'nap', '.gt.')
      endif
      if (switch%b2srdt_phm0.ne.1.0_R8) then
        write(*,*)
     & 'Warning: b2srdt_phm0 should be 1.0 for time-dependent mode'
      endif

* ..include the contributions from dtim
*   ..modify snadt
      do iCv = 1, nCi !only internal cells
        ttim=dtim*dtco(isb,mpg%cvReg(iCv))*time_factor(iCv)
        t0=(switch%b2srdt_phm0/ttim)*cvVol(iCv)                                     !srv 06.01.21
        snadt(iCv,0) = t0*nap(iCv)                                  !srv 06.01.21
        snadt(iCv,1) =-t0                                           !srv 06.01.21
      enddo


!WG_TODO cdjm Store the total change in linearised sources for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO        b2srdt_sna0to1(:,:,:,isb) = snadt
!WG_TODO      endif

*    The diagnostic output of the time derivative will be done
*    in the subroutines where the solution is calculated.
      if (switch%b2scdt_iout.eq.1) then
        write (chns,'(i3.3)') isb
        do k=0,1
         write (chk,'(i1)') k
         call my_out_us(70,nCv,0,snadt(:,k),
     &                       'b2scdt_snadt'//chk//chns)
        enddo
      endif                                                           !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2scdt

      end subroutine b2scdt
