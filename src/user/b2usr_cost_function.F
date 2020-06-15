      subroutine b2usr_cost_function(mpg, state, J)
      use b2mod_types
      use b2mod_constants
      use b2mod_user_namelist
      use b2us_map
      use b2us_plasma
      implicit none
      type(mapping), intent (in)    :: mpg
      type(B2State), intent (inout)    :: state
      real(kind=R8) J
      integer icf, ic1, ic2, iCv

      J = 0.0_R8

      if (ncf > 0) then
        do icf = 1, ncf
          ic1 = mpg%cfCvP(icf,1)
          ic2 = ic1 + mpg%cfCvP(icf,2) - 1
          do iCv = ic1, ic2
            J = max(J,state%pl%te(mpg%cfCv(iCv)))
          enddo
        enddo
      J= J/ev
      endif     

      return
      end subroutine b2usr_cost_function
