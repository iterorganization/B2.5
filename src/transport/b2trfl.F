*.specification

      subroutine b2trfl (iCv, switch, mpg, Ainv, zeffne, amti, zetap,
     &  gamma, cF_lim, cF_limh, cF_limB, cF_limBh)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_constants
      use b2us_map
      use b2mod_switches
      use b2mod_zhfrtf
     , , only: gtalc, nal
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer iCv
      real (kind=R8) :: gamma                                           !som 20.03.21
      real (kind=R8) :: zeffne, zetap, amti(0:nspecies-1),
     &  Ainv(0:2*nspecies-1,0:2*nspecies-1)
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  cF_lim(0:nspecies-1),
     &  cF_limh(0:nspecies-1),
     &  cF_limB(0:nspecies-1),
     &  cF_limBh(0:nspecies-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRFL computes the flux-limiting coefficients.
*
*  3. description
*
*     This routine computes:
*     - cF_lim: the flux-limiting coefficient for the thermal force coefficient
*     - cF_limh: the flux-limiting coefficient for the heat flux
*     - cF_limB: the flux-limiting coefficient for the charge state correction
*     - cF_limBh: the flux-limiting coefficient for the charge state correction for the heat flux
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer i
      real (kind=R8) :: coef,t0,t1,t2
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trfl')
*   ..main
      if (cflim(4).ne.0.0_R8.and..not.mpg%cvOnClosedSurface(iCv)) then
        coef=2.5_R8*zetap/mp/(zeffne**2)
        do i=0, nspecies-1
          t0=abs(sum(Ainv(i,0:nspecies-1)*nal(iCv,0:nspecies-1)*
     &     gtalc(iCv,0:nspecies-1)))
          t1=switch%zhflcorr*cflim(4)/sqrt(amti(i))
          t2=switch%zhflcorrh*switch%zhflcorr*cflim(4)/sqrt(amti(i))

          cF_lim(i)=
     &     1.0_R8/(1.0_R8+(coef*t0/t1)**gamma)**(1.0_R8/gamma)
          cF_limh(i)=
     &     1.0_R8/(1.0_R8+(coef*t0/t2)**gamma)**(1.0_R8/gamma)
          cF_limB(i)=0.5_R8+0.5_R8*
     &     tanh(10.0_R8*(t1/(coef*t0)-1.0_R8))
          cF_limBh(i)=0.5_R8+0.5_R8*
     &     tanh(10.0_R8*(t2/(coef*t0)-1.0_R8))
        enddo
      else
        cF_lim(0:nspecies-1)=1.0_R8
        cF_limh(0:nspecies-1)=1.0_R8
        cF_limB(0:nspecies-1)=1.0_R8
        cF_limBh(0:nspecies-1)=1.0_R8
      endif

      call subend ()
      return

      end subroutine b2trfl
