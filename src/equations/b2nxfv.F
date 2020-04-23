*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfv (nCv, nFc, nVx, isb, geo, mpg, cvsbhz_cl, flubv)            !srv 19.01.05 01.07.05 11.09.09 { 16.10.17
      use b2mod_types
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   cvsbhz_cl(nFc,0:1)               !srv 16.10.17
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   flubv(nFc,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXFV calculates the convective term due to new form of ion  
*     viscosity term in the momentum correction  equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall
      real (kind=R8) ::  phm0 = 1.0_R8                                  !srv 06.12.18
      real (kind=R8) ::
     *   ft, dloghz(nFc)                                                !srv 16.10.17
      save ncall
      external subini, subend, xertst, ipgeti, ipgetr, b2scopy
      external b2xvff, b2xvfx, b2xvfv, b2xvsg
*   ..initialization
      data ncall/0/

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
      if (ncall.eq.0) then                                              !srv 06.12.18 {
       call ipgetr ('b2nxfv_phm0', phm0)
      endif                                                             !srv 06.12.18 }
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
      endif

* ..compute convection term
      ft = 4.0_R8/3.0_R8
*   ..compute flubvx
      if(nint(zn(isb)).eq.1 .and. nint(zamax(isb)).eq.1) then         !srv 26.10.08
        call diff_p(nCv,nFc,nVx,1,geo,mpg,log(geo%cvHz),log(geo%vxHz),
     .    dloghz)
        flubv(:,0) = phm0*ft*cvsbhz_cl(:,0)*dloghz                    !srv 06.12.18
        flubv(:,1) = 0.0_R8
      else
        flubv = 0.0_R8
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxfv

      end                                                            !srv 19.01.05 01.07.05 11.09.09 }
