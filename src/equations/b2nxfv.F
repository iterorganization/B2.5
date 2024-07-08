*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfv (nCv, nFc, nVx, isb, ismain, switch, geo, mpg,  !srv 19.01.05 01.07.05 11.09.09 { 16.10.17
     &  rho_rel, cvsbhz_cl, cvsbhz, flubv)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8), intent(in) ::
     *   rho_rel(nCv), cvsbhz_cl(nFc,0:1), cvsbhz(nFc,0:1)              !srv 16.10.17
*   ..output arguments (unspecified on entry)
      real (kind=R8), intent(out) ::
     *   flubv(nFc,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXFV calculates the convective term due to new form of ion  
*     viscosity term in the momentum correction equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iVx
      real (kind=R8) ::
     *   ft, dloghz(nFc,0:1), dlogrho(nFc), wrk(nCv), wrkv(nVx)         !srv 16.10.17
      external xertst
      external b2xvff, b2xvfx, b2xvfv, b2xvsg
*   ..initialization

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nxfv')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')

* ..compute convection term
      ft = 4.0_R8/3.0_R8
*   ..compute flubvx
      if (isb.eq.ismain.or.
     &    (switch%Zhdanov_closure.ne.0 .and.
     &     switch%Zhdanov_test.ne.1 .and.
     &     switch%Zhdanov_visu.ne.0 .and.
     &     .not.is_neutral(isb))) then                                  !srv 26.10.08
        ! compute d(logHz)
        wrk = log(geo%cvHz)
        do iVx = 1, nVx
          if (geo%vxHz(iVx).gt.0.0_R8) then
            wrkv(iVx) = log(geo%vxHz(iVx))
          else
            wrkv(iVx) = 0.0_R8
          end if
        end do
        call diff_p(nCv,nFc,nVx,1,geo,mpg,wrk,wrkv,dloghz(1,0))
        ! compute d(logRho)
        if (switch%Zhdanov_closure.ne.0 .and.
     &      switch%Zhdanov_test.ne.1 .and.
     &      switch%Zhdanov_visu.ne.0 .and.
     &      switch%zh_vis_type.eq.1) then
          wrk = log(rho_rel)
          call diff_p(nCv,nFc,nVx,0,geo,mpg,wrk,wrkv,dlogrho)
        else
          dlogrho = 0.0_R8
        endif
        flubv(:,0) = switch%b2nxfv_phm0*ft*cvsbhz_cl(:,0)*
     &    (dloghz(:,0) + dlogrho)  !srv 06.12.18
        flubv(:,1) = 0.0_R8
      elseif (nint(zn(isb)).eq.1 .and. nint(zamax(isb)).eq.0) then
        wrk = log(geo%cvHz)
        do iVx = 1, nVx
          if (geo%vxHz(iVx).gt.0.0_R8) then
            wrkv(iVx) = log(geo%vxHz(iVx))
          else
            wrkv(iVx) = 0.0_R8
          end if
        end do
        call diff(nCv,nFc,nVx,1,geo,mpg,wrk,wrkv,dloghz)
        flubv(:,0) = switch%b2nxfv_phm0*switch%b2nxfv_phm1
     &                                      *cvsbhz(:,0)*dloghz(:,0)
        flubv(:,1) = switch%b2nxfv_phm0*switch%b2nxfv_phm1
     &                                      *cvsbhz(:,1)*dloghz(:,1)
      else
        flubv = 0.0_R8
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxfv

      end subroutine b2nxfv                                             !srv 19.01.05 01.07.05 11.09.09 }
