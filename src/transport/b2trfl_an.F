*.specification

      subroutine b2trfl_an (ix, iy, c_hTA_an, amnucl, amti, zetap,
     &  gamma, cF_lim, cF_limh, cF_limB, cF_limBh)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_indirect
      use b2mod_constants
      use b2mod_zhfrtf
     , , only: gtalc, nal, z2n_xy
      use b2mod_b2zhco
     , , only: zhflcorr, zhflcorrh
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ix,iy
      real (kind=R8) :: gamma                                           !som 20.03.21
      real (kind=R8) :: amnucl(0:nspecies-1), zetap,
     &  amti(0:nspecies-1), c_hTA_an(0:nspecies-1)
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
*     B2TRFL_AN computes the flux-limiting coefficients using
*     analytically calculated ion heat conductivity.
*
*  3. description
*
*     This routine computes:
*     - cF_lim: the flux-limiting coefficient for the thermal force coefficient
*     - cF_limh: the flux-limiting coefficient for the heat flux
*     - cF_limB: the flux-limiting coefficient for the charge state correction
*     - cF_limBh: the flux-limiting coefficient for the charge state correction for the heat flux
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer i
      real (kind=R8) :: coef,t0,t1,t2
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trfl_an')
*   ..main
      if (cflim(4).ne.0.0_R8.and..not.on_closed_surface(ix,iy)) then
        coef=zetap/mp
        do i=0, nspecies-1
          t0=nal(ix,iy,i)/z2n_xy(ix,iy,i)**2*sqrt(2.0_R8/amnucl(i))*
     &     abs(c_hTA_an(i)*gtalc(ix,iy,i))
          t1=zhflcorr*cflim(4)/sqrt(amti(i))
          t2=zhflcorrh*zhflcorr*cflim(4)/sqrt(amti(i))

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

      end subroutine b2trfl_an
