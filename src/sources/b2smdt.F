*-----------------------------------------------------------------------
*.specification

      subroutine b2smdt (nCv, nCi, ns, dtim, hz,                     !srv 09.01.01 01.07.08
     &  cvVol, na, uap, switch, mpg, smodt, isb)                     !srv 11.09.09
      use b2mod_types
      use b2mod_numerics_namelist
      use b2mod_constants
      use b2mod_b2cmpa
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     1 , only : b2srdt_smo0to3, balance_netcdf
      use b2mod_switches
      use b2us_map
      use b2mod_subsys
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2smdt
      implicit none
*   ..input arguments (unchanged on exit)
      type (switches), intent (in) :: switch
      type (mapping), intent (in) :: mpg
      integer nCv, nCi, ns, isb
      real (kind=R8), intent(in) ::
     *  dtim, cvVol(nCv), na(nCv,0:ns-1),
     &  uap(nCv,0:ns-1), hz(nCv)                                     !srv 09.01.01
*   ..input/output arguments
      real (kind=R8) ::
     *  smodt(nCv,0:3,0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SMDT modifies the source coefficients to account for the
*     implicit time discretisation for velocity equations.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, k
      character chns*3, chk*1
      real (kind=R8) :: t0, ttim
*   ..procedures
      external xertst, ipgetr
      external b2xvsg

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2smdt')
*   ..test sign of dtim
      call xertst (0.0_R8.lt.dtim, 'faulty argument dtim')
*   ..extensive tests on first few calls
      if (ncall_b2smdt.lt.3) then
*    ..test old state
       call b2xvsg (nCv*ns, na, 1, 'na', '.gt.')
      endif
      if (switch%b2srdt_phm1.ne.1.0_R8) then
        write(*,*)
     & 'Warning: b2srdt_phm1 should be 1.0 for time-dependent mode'
      endif

*   ..modify time-stepping sources according to time integration scheme
*   .. ts_factor = time-stepping multiplication factor
* ..include the contributions from dtim

*   ..modify smodt
      do iCv = 1, nCi !only internal cells
        ttim=dtim*dtmo(isb,mpg%cvReg(iCv))*time_factor(iCv)
        t0 = (switch%b2srdt_phm1/ttim)*cvVol(iCv)*hz(iCv)
     &       *ts_factor
        smodt(iCv,0,isb) = t0*am(isb)*mp*uap(iCv,isb)*na(iCv,isb)        !srv 06.01.21 {
        smodt(iCv,1,isb) = 0.0_R8
        smodt(iCv,2,isb) = 0.0_R8
        smodt(iCv,3,isb) = -t0                                           !srv 06.01.21 }
      enddo

cdjm Store the total change in linearised sources for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO         b2srdt_smo0to3 = smodt
!WG_TODO       endif

*    The diagnostic output of the time derivative will be done
*    in the subroutines where the solution is calculated.
      if (switch%b2smdt_iout.eq.1) then
        write (chns,'(i3.3)') isb
        do k=0,3
         write (chk,'(i1)') k
         call my_out_us(70,nCv,0,smodt(:,k,isb),
     &                       'b2smdt_smodt'//chk//chns)
        enddo
      endif                                                          !srv 11.09.09 }

* ..return
      ncall_b2smdt = ncall_b2smdt+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2smdt

      end subroutine b2smdt
