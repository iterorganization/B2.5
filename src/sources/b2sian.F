*-----------------------------------------------------------------------
*.specification

      subroutine b2sian (nCv, nFc, isb, geo, mpg, fchanml, smban)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     & , only : b2sian_smo0to3, balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, isb
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
*     to the anomalous current
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv, k
      integer, save :: iout = 0
      character chns*3, chark*1
      real (kind=R8) ::
     *    phm0, eps, wrk(nCv), wrkf(nFc)
      save ncall, phm0, eps
*   ..procedures
      external subini, subend, xertst, ipgeti, ipgetr
      external b2xvff, b2xvsg
*   ..initialization
      data ncall/0/, phm0/1.0e0_R8/, eps/1.0e-60_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sian')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2sian_phm0', phm0)
       call ipgeti ('b2sian_iout', iout)
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      smban = 0.0_R8
      if(nint(zn(isb)).eq.1.and.nint(zamax(isb)).eq.1) then
*   ..  compute smban
        wrkf = fchanml(:,1)/(geo%fcS*geo%fcQalf(:,1) + eps)
        call intcell(nFc,nCv,mpg,mpg%intcellR,wrkf,wrk)
        do iCv = 1, mpg%nCi
            smban(iCv,0) = -phm0*geo%cvH(iCv)*geo%cvVol(iCv)*
     &          geo%cvBb(iCv,0)/geo%cvBb(iCv,2)*geo%cvBb(iCv,3)*wrk(iCv)
        enddo
      endif
cdjm Store linearised source for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO        b2sian_smo0to3(-1:nx,-1:ny,0:3,isb)=smban
!WG_TODO      endif

      if (iout.ne.0) then
        write (chns,'(i3.3)') isb
        do k=0,3
          write (chark,'(i1)') k
          call my_out_us(70,nCv,smban(1,k),
     .                         'b2sian_smoan'//chark//chns)
        enddo
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sian

      end subroutine b2sian
