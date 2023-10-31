*-----------------------------------------------------------------------
*.specification

      subroutine b2smdt (nCv, nCi, ns, dtim, hz,          !srv 09.01.01 01.07.08
     &  cvVol, na0, ua0, switch, mpg, smodt)              !srv 11.09.09
      use b2mod_types
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_balance !djm Jan2017
     1 , only : b2srdt_smo0to3, balance_netcdf
      use b2mod_switches
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      type (switches), intent (in) :: switch
      type (mapping), intent (in) :: mpg
      integer nCv, nCi, ns
      real (kind=R8), intent(in) ::
     *  dtim, cvVol(nCv), na0(nCv,0:ns-1),
     &  ua0(nCv,0:ns-1), hz(nCv)                                     !srv 09.01.01
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
      integer ncall, iCv, is, k
      character chns*3, chk*1
      real (kind=R8) :: t0, ttim
*   ..procedures
      external xertst, ipgetr
      external b2xvsg
*   ..initialisation
      data ncall/0/                                       !xpb

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2smdt')

*   ..test sign of dtim
      call xertst (0.0_R8.lt.dtim, 'faulty argument dtim')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test old state
       call b2xvsg (nCv*ns, na0, 1, 'na0', '.gt.')
      endif

* ..include the contributions from dtim

*   ..modify smodt
      do is = 0, ns-1
        do iCv = 1, nCi !only internal cells
            ttim=dtim*dtmo(is,mpg%cvReg(iCv))*time_factor(iCv)
            t0 = (switch%b2srdt_phm1/ttim)*cvVol(iCv)*hz(iCv)                        !srv 09.01.01
            smodt(iCv,0,is) = t0*am(is)*mp*ua0(iCv,is)*na0(iCv,is) !srv 06.01.21 {
            smodt(iCv,1,is) = 0.0_R8
            smodt(iCv,2,is) = 0.0_R8
            smodt(iCv,3,is) = -t0                                      !srv 06.01.21 }
        enddo
       enddo

cdjm Store the total change in linearised sources for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO         b2srdt_smo0to3 = smodt
!WG_TODO       endif

*    The diagnostic output of the time derivative will be done
*    in the subroutines where the solution is calculated.
      if (switch%b2smdt_iout.eq.1) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        do k=0,3
         write (chk,'(i1)') k
         call my_out_us(70,nCv,smodt(:,k,is),
     &                       'b2smdt_smodt'//chk//chns)
        enddo
       enddo
      endif                                                           !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2smdt

      end subroutine b2smdt
