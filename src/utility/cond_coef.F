      function cond_coef(inertial,T)
      use b2mod_types
      implicit none
      real(kind=R8) :: cond_coef
      logical inertial
      real(kind=R8) :: T, Tc

      Tc = T - 273.15_R8

      if (inertial) then
! we assume a residual contact heat conductivity of 100 W.m-2.K-1
        cond_coef = 100.0_R8
      else
! The following are fits to the result obtained by Nukiyama for the heat conduction
! transfer coefficient between the CuCrZr pipe and the water coolant,
! applicable to the ITER divertor targets
! Full reference can be found at: https://user.iter.org/?uid=35L29Y
!      cond_coef = 1.64908e-8_R8*Tc**6 - 1.24649e-5_R8*Tc**5
!     @       + 3.4938e-3_R8*Tc**4 - 0.44568_R8*Tc**3 + 25.5532_R8*Tc**2
!     @       - 502.333_R8*Tc + 77703.2_R8
!      cond_coef = 1.357954e-3_R8*Tc**3 + 6.616505e-1_R8*Tc**2 - 24.41418_R8*Tc
!     @       + 7.725710e4_R8
        cond_coef = 8.0e4_R8
      end if

      return
      end function cond_coef
