      module b2mod_transport_fun                                        !srv 13.01.17
      use b2mod_types

      contains

      function fce1(z)
      implicit none
      real (kind=R8) :: z, fce1
      fce1=((1.0_R8+0.24_R8*z)*(1.0_R8+0.93_R8*z))/
     &     ((1.0_R8+2.56_R8*z)*(1.0_R8+0.29_R8*z))
      end function fce1

      function fce2(z)
      implicit none
      real (kind=R8) :: z, fce2

      fce2=((1.0_R8+1.40_R8*z)*(1.0_R8+0.52_R8*z))/
     &     ((1.0_R8+2.56_R8*z)*(1.0_R8+0.29_R8*z))*1.56_R8
      end function fce2

      function fce2n(z)
      implicit none
      real (kind=R8) :: z, fce2n

      fce2n=fce2(z)/(z+sqrt(2.0_R8)/2.0_R8)
      end function fce2n

      function fse_ce1(z)
      implicit none
      real (kind=R8) :: z, fse_ce1
      fse_ce1=fce1(z)
      end function fse_ce1

      function fal_cen(z)
      implicit none
      real (kind=R8) :: z, fal_cen
      fal_cen = -fce2n(z)/fce1(z)
      end function fal_cen

      end module b2mod_transport_fun                                     !srv 13.01.17
