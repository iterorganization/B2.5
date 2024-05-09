*-----------------------------------------------------------------------
*.specification

      subroutine b2sian(nCv, nFc, isb, switch, geo, mpg, fchanml, smban)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
#ifndef FIXED_POINT
      use b2mod_balance !djm Jan2017
     & , only : b2sian_smo0to3, balance_netcdf
#endif
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : b2sian_eps
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, isb
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &   fchanml(nFc,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smban(nCv,0:3)
*   ..workspace
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIAN computes contribution to the parallel momentum source due
*     to the anomalous current.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, k
      character chns*3, chark*1
      real (kind=R8) ::
     *    wrk(nCv), wrkf(nFc)
*   ..procedures
      external xertst
      external b2xvff, b2xvsg
*   ..initialization

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sian')
*   ..set internal parameters on first call
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      smban = 0.0_R8
      if(nint(zn(isb)).eq.1.and.nint(zamax(isb)).eq.1) then
*   ..  compute smban
        wrkf = fchanml(:,1)/(geo%fcS*geo%fcQalf(:,1) + b2sian_eps)
        call intcell(nFc,nCv,mpg,mpg%intcellR,wrkf,wrk)
        do iCv = 1, mpg%nCi
            smban(iCv,0) = -switch%b2sian_phm0*
     &          geo%cvHz(iCv)*geo%cvVol(iCv)*
     &          geo%cvBb(iCv,0)/geo%cvBb(iCv,2)*geo%cvBb(iCv,3)*wrk(iCv)
        enddo
      endif
#ifndef FIXED_POINT
cdjm Store linearised source for balance
      if (balance_netcdf.ne.0) then
        b2sian_smo0to3(1:nCv,0:3,isb)=smban
      endif
#endif

      if (switch%b2sian_iout.ne.0) then
        write (chns,'(i3.3)') isb
        do k=0,3
          write (chark,'(i1)') k
          call my_out_us(70,nCv,0,smban(1,k),
     .                         'b2sian_smoan'//chark//chns)
        enddo
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sian

      end subroutine b2sian
