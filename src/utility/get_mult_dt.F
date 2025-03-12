*-----------------------------------------------------------------------
*.specification

      subroutine get_mult_dt(style, nCv, te, te_min, te_max, 
     &                       mult_max, mult_dt)
      use b2mod_types
      use b2mod_constants
      implicit none
*   ..input arguments
      integer, intent(in) :: style, nCv
      real (kind=R8), intent(in) ::
     *   te_min, te_max, mult_max, te(nCv)
*   ..output arguments
      real (kind=R8), intent(out) ::
     *   mult_dt(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*     get_mult_dt calculates a local (cell-based) multiplier to the timestep. 
*     The style switch determines the model.
*
*     style = 0: no multiplier, dt_mult = 1.0
*
*     style = 1: based on Te. In low Te regions (typically more interaction with neutrals) a small timestep is preferred. 
*                In high Te regions there is typically less neutral interaction, or this interaction
*                is a least less critical from a numerical point of view. 
*
*     style = 2: based on the magnitude of the sources compared to the particle/momentum/energy content of the cell. 
*
*     Other styles can be added in the future.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: iCv
      real (kind=R8) :: te_min_J, te_max_J, exp1
*   ..procedures
      external xertst
*-----------------------------------------------------------------------
*.computation  


      if (style.eq.0) then
        mult_dt = 1.0_R8
      elseif (style.eq.1) then
        te_min_J = te_min*ev !! convert from eV to joule
        te_max_J = te_max*ev
        do iCv = 1, nCv
          if (te(iCv).ge.te_max_J) then
            mult_dt(iCv) = mult_max
          elseif (te(iCv).le.te_min_J) then
            mult_dt(iCv) = 1.0_R8
          else
            exp1 = log10(mult_max)*(log10(te(iCv))-log10(te_min_J)) / 
     &                              (log10(te_max_J)-log10(te_min_J))
            mult_dt(iCv) = 10**exp1
          endif
        enddo
      elseif (style.eq.2) then !! This option is worked out in b2srdt itself.
        mult_dt = 1.0_R8 
      else
        call xerrab('This b2mndt_var_style is not available')
      endif

*   ..return
      return

*-----------------------------------------------------------------------


      end subroutine get_mult_dt
