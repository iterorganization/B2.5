*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfv (nCv, nFc, nVx, isb, switch, geo, mpg,        !srv 19.01.05 01.07.05 11.09.09 { 16.10.17
     &  cvsbhz_cl, cvsbhz, flubv)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8), intent(in) ::
     *   cvsbhz_cl(nFc,0:1), cvsbhz(nFc,0:1)                            !srv 16.10.17
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
      real (kind=R8) ::
     *   ft, dloghz(nFc,0:1), logCvHz(nCv), logVxHz(nVx)                !srv 16.10.17
      external xertst, b2scopy
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
      if (nint(zn(isb)).eq.1 .and. nint(zamax(isb)).eq.1) then          !srv 26.10.08
        logCvHz = log(geo%cvHz)
        logVxHz = log(geo%vxHz)
        call diff_p(nCv,nFc,nVx,1,geo,mpg,logCvHz,logVxHz,dloghz(1,0))
        flubv(:,0) = switch%b2nxfv_phm0*ft*cvsbhz_cl(:,0)*dloghz(:,0)                 !srv 06.12.18
        flubv(:,1) = 0.0_R8
      elseif (nint(zn(isb)).eq.1 .and. nint(zamax(isb)).eq.0) then
        logCvHz = log(geo%cvHz)
        logVxHz = log(geo%vxHz)
        call diff(nCv,nFc,nVx,1,geo,mpg,logCvHz,logVxHz,dloghz)
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

      end subroutine b2nxfv                                           !srv 19.01.05 01.07.05 11.09.09 }
