*.specification

      subroutine b2trfl_vis (iCv, mpg, ns, eta, pa, gamma, cfl, cF_limv)
      use b2mod_types
      use b2mod_zhfrtf , only: gavm_uc
      use b2mod_subsys
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer iCv,ns
      real (kind=R8) :: gamma, cfl                                      !som 20.03.21
      real (kind=R8) :: eta(0:ns-1), pa(0:ns-1)
      type (mapping), intent (in) :: mpg
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  cF_limv(0:ns-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRFL_VIS computes the flux-limiting coefficient for the stress
*     viscosity.
*
*  3. description
*
*     This routine computes:
*     - cF_limv: the flux-limiting coefficient for the thermal force
*                coefficient
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer is
      real (kind=R8) :: t0,t1
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trfl_vis')
*   ..main
      if (cfl.ne.0.0_R8.and..not.mpg%cvOnClosedSurface(iCv)) then
        do is=0, ns-1
          t0=abs(gavm_uc(iCv))
          t1=cfl*pa(is)

          cF_limv(is)=
     &     1.0_R8/(1.0_R8+(eta(is)*t0/t1)**gamma)**(1.0_R8/gamma)
        enddo
      else
        cF_limv(0:ns-1)=1.0_R8
      endif

      call subend ()
      return

      end subroutine b2trfl_vis
