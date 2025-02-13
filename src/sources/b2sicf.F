*-----------------------------------------------------------------------
*.specification

      subroutine b2sicf (nCv, nFc, nVx, isb, mb, switch, geo, mpg,     !srv 09.07.01 28.01.02 22.07.05
     *                  nb, ub, smbcf, ctcf, wrk0)                     !srv 06.04.07
      use b2mod_types
      use b2mod_constants
      use b2mod_switches
      use b2us_geo
      use b2us_map
#ifndef FIXED_POINT
      use b2mod_balance !djm Jan2017
     & , only : b2sicf_smo0to3, balance_netcdf
#endif
      use b2mod_ad
     & , only : ncall_b2sicf
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb                                        !srv 28.01.02
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   mb, nb(nCv), ub(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smbcf(nCv,0:3), ctcf(nCv)                                      !srv 06.04.07
*   ..workspace
      real (kind=R8) ::
     *   wrk0(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SICF computes centrifugal force
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, k
      character chns*3, chark*1                                         !srv 17.06.02
      real (kind=R8) :: wrks(nCv)
*   ..procedures
      external xertst
      external b2xvsg

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sicf')
*   ..set internal parameters on first call
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test mb
      call xertst (0.le.mb, 'faulty argument mb')
*   ..extensive tests on first few calls
      if (ncall_b2sicf.lt.3) then
*    ..test sign of nb
       call b2xvsg (nCv, nb, 1, 'nb', '.gt.')
      endif

* ..compute centrifugal force
      smbcf = 0.0_R8
      ctcf = 0.0_R8                                                     !srv 06.04.07
      if (switch%b2sicf_style.eq.0) then
        call gradc_p(nCv,nFc,nVx,1,geo,mpg,geo%cvHz,geo%vxHz,wrk0)

*   .. compute smbcf
        do iCv = 1, mpg%nCi
          smbcf(iCv,0) = switch%b2sicf_phm0*geo%cvVol(iCv)*             !srv factor hz/hz=1
     &     (geo%cvBb(iCv,0)/geo%cvBb(iCv,3))*mb*ub(iCv)**2              !srv 05.04.07
     &      *nb(iCv)*wrk0(iCv)
          ctcf(iCv) = switch%b2sicf_phm0*switch%b2sicf_phm1*
     &                                        geo%cvVol(iCv)*           !srv 06.04.07 {
     &     (geo%cvBb(iCv,0)/geo%cvBb(iCv,3))*mb*2.0_R8*ub(iCv)
     &      *nb(iCv)*wrk0(iCv)                                          !srv 06.04.07 }
        enddo

      else if (switch%b2sicf_style.eq.1) then
        call xerrab ('b2sicf -- style 1 not for WG')
      endif

#ifndef FIXED_POINT
cdjm Store linearised source for balance
      if (balance_netcdf.ne.0) then
        b2sicf_smo0to3(1:nCv,0:3,isb)=smbcf
      endif
#endif

      if (switch%b2sicf_iout.ne.0) then                                 !srv 17.06.02 {
        write (chns,'(i3.3)') isb
        do k=0,3
          write (chark,'(i1)') k
          call my_out_us(70,nCv,0,smbcf(1,k),
     .                         'b2sicf_smocf'//chark//chns)
        enddo
        wrks=smbcf(:,0)+smbcf(:,1)*ub+
     &       smbcf(:,2)*nb*mb+smbcf(:,3)*ub*nb*mb
        call my_out_us(70,nCv,0,wrks,'b2sicf_smocf'//chns)
        call my_out_us(70,nCv,0,ctcf,'b2sicf_ctcf'//chns)
      endif                                                             !srv 17.06.02 }

* ..return
      ncall_b2sicf = ncall_b2sicf+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sicf

      end subroutine b2sicf
