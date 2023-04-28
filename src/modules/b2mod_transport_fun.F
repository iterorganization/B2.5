      module b2mod_transport_fun                                        !srv 13.01.17
      use b2mod_types

      real (kind=R8) c11e, c13e, c15e
      parameter (c11e = 1.0000_R8, c13e = 0.9487_R8, c15e =-0.8964_R8)
      real (kind=R8) c31e
      parameter (c31e = c13e)

      contains

      function c33e(z)
      implicit none
      real(kind=R8) :: c33e, z
      c33e = 1.3000_R8 + 0.5657_R8/z
      return
      end function c33e

      function c35e(z)
      implicit none
      real(kind=R8) :: c35e, z
      c35e = -1.3040_R8 - 0.3207_R8/z
      return
      end function c35e

      function c53e(z)
      implicit none
      real(kind=R8) :: c53e, z
      c53e = c35e(z)
      return
      end function c53e

      function c55e(z)
      implicit none
      real(kind=R8) :: c55e, z
      c55e = 1.5464_R8 + 0.9091_R8/z
      return
      end function c55e

      function d13e(z)
      implicit none
      real(kind=R8) :: d13e, z
      d13e = c11e   *c33e(z) - c13e   *c13e
      return
      end function d13e

      function d15e(z)
      implicit none
      real(kind=R8) :: d15e, z
      d15e = c11e   *c55e(z) - c15e   *c15e
      return
      end function d15e

      function d35e(z)
      implicit none
      real(kind=R8) :: d35e, z
      d35e = c33e(z)*c55e(z) - c35e(z)*c35e(z)
      return
      end function d35e

      function f135e(z)
      implicit none
      real(kind=R8) :: f135e, z
      f135e = c55e(z)*d13e(z) - c35e(z)*h315e(z) - c15e   *h135e(z)
      return
      end function f135e

      function fal(z)
      implicit none
      real(kind=R8) :: fal, z
      fal = -1.25_R8*z/(1.0_R8+0.41_R8*z)
      return
      end function fal

      function fal_cen(z)
      implicit none
      real (kind=R8) :: z, fal_cen
      fal_cen = -fce2n(z)/fce1(z)
      end function fal_cen

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

      function fke_Zh(z)
      implicit none
      real(kind=R8) :: fke_Zh, z
      fke_Zh = (3.9_R8+2.3_R8/z)/(0.31_R8+1.2_R8/z+0.41_R8/z**2)
      return
      end function fke_Zh

      function fke(z)
      implicit none
      real(kind=R8) :: fke, z
      fke = 2.16_R8*z/(1.0_R8+0.27_R8*z)
      return
      end function fke

      function fse(z)
      implicit none
      real(kind=R8) :: fse, z
      fse = 4.95_R8*z/(1.0_R8+1.54_R8*z)
      return
      end function fse

      function fse_ce1(z)
      implicit none
      real (kind=R8) :: z, fse_ce1
      fse_ce1=fce1(z)
      end function fse_ce1

      function h135e(z)
      implicit none
      real(kind=R8) :: h135e, z
      h135e = c15e   *c33e(z) - c13e   *c35e(z)
      return
      end function h135e

      function h153e(z)
      implicit none
      real(kind=R8) :: h153e, z
      h153e = c13e   *c55e(z) - c15e   *c53e(z)
      return
      end function h153e

      function h315e(z)
      implicit none
      real(kind=R8) :: h315e, z
      h315e = c35e(z)*c11e    - c31e   *c15e
      return
      end function h315e

      end module b2mod_transport_fun                                     !srv 13.01.17
