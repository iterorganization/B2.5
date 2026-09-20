*-----------------------------------------------------------------------
*.specification

      subroutine b2sicf (nCv, nFc, nVx, isb, mb, switch, geo, mpg,      !srv 09.07.01 28.01.02 22.07.05
     *                   nb, ub, smbcf, ctcf, wrk0)                     !srv 06.04.07
      use b2mod_types
      use b2mod_constants
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_openmp
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
*     ..the centrifugal force is quadratic in the (unknown) parallel
*     ..velocity; Newton-linearise ub**2 around the current ub instead
*     ..of dumping the full quadratic into the constant slot (0) with
*     ..no matrix feedback (slots 1/3 left at their zero fill) -- the
*     ..latter left the term fully explicit and able to run away once
*     ..ub got large, with no compensating diagonal contribution.
*     ..ub**2 = ub_old*(2*ub - ub_old) + O((ub-ub_old)**2), so put the
*     ..-ub_old**2 piece in the rob-weighted constant slot (2) and the
*     ..+2*ub_old*ua piece in the rob*ua-weighted slot (3); mb*nb=rob
*     ..is supplied by the caller's reconstruction, so it must not be
*     ..included again here. The linear (slot-3) piece is only used
*     ..when it is stabilising (<=0), matching the codebase invariant
*     ..that smb(:,:,3) must never be positive -- this mirrors what the
*     ..old one-sided ctcf=min(0,ctcfb) diagonal patch did, just routed
*     ..through the smbcf(2:3) slots. When the derivative would
*     ..destabilise (>0), fall back to the fully explicit treatment.
        do iCv = 1, mpg%nCi
          wrks(iCv) = switch%b2sicf_phm0*geo%cvVol(iCv)*                !srv factor hz/hz=1
     &     (geo%cvBb(iCv,0)/geo%cvBb(iCv,3))*wrk0(iCv)                  !srv 05.04.07
          if (wrks(iCv)*ub(iCv).le.0.0_R8) then
            smbcf(iCv,2) = -wrks(iCv)*ub(iCv)**2
            smbcf(iCv,3) =  wrks(iCv)*2.0_R8*ub(iCv)
          else
            smbcf(iCv,2) =  wrks(iCv)*ub(iCv)**2
            smbcf(iCv,3) =  0.0_R8
          endif
c         ..ctcf left at zero (see init above): now handled through
c         ..smbcf(2:3) above, not the one-sided ctcf patch
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
#ifdef _OPENMP
!$OMP MASTER
        write(*,*) 'B2SICF OpenMP warning: no file output in ',
     &             'parallel mode'
!$OMP END MASTER
#else
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
#endif
      endif                                                             !srv 17.06.02 }

* ..return
      ncall_b2sicf = ncall_b2sicf+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sicf

      end subroutine b2sicf
