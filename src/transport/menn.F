
      subroutine menn(ta,tb,ma,mb,mm,nn)
c--------------------------------------------------------------------
c     subroutine that calculates the friction coefficients m^ab and 
c     n^ab 
c
c     input    ta   : the temperature of species a
c              tb   : the temperature of species b
c              ma   : the mass of species a
c              mb   : the mass of species b
c     output   mm   : array(3,3) containing the test particle friction
c                     coefficients
c              nn   : array(3,3) containing the field particle friction
c                     coefficients
c--------------------------------------------------------------------

      use b2mod_types

      implicit none

      real (kind=R8) ::
     *        ta, tb, ma, mb, mm(3,3), nn(3,3), xab, tab, mab

      tab = ta / tb
      mab = ma / mb
      xab = mab / tab

      mm(1,1) = - (1.0_R8 + mab) / (1. + xab)**1.5_R8
      mm(1,2) =   1.5_R8*(1 + mab)/(1. + xab)**2.5_R8
      mm(2,1) =   mm(1,2)
      mm(2,2) = - (3.25_R8 + 4.0_R8*xab + 7.5_R8*xab**2) /
     .            (1.0_R8 + xab)**2.5_R8
      mm(1,3) = - 15.0_R8 *(1.0_R8+mab)/(8.0_R8*(1.+xab)**3.5_R8)
      mm(3,1) =   mm(1,3)
      mm(2,3) =   (69.0_R8/16.0_R8+6.0_R8*xab + 15.75_R8*xab**2)/
     .            (1.0_R8+xab)**3.5_R8
      mm(3,2) =   mm(2,3)
      mm(3,3) = - (433.0_R8/64.0_R8+ 17.0_R8*xab +
     .             459.0_R8*(xab**2)/8.0_R8 + 28.0_R8*xab**3 +
     +             175.0_R8*(xab**4)/8.0_R8)/(1.0_R8 + xab)**4.5_R8

      nn(1,1) = - mm(1,1)
      nn(2,1) = - mm(2,1)
      nn(3,1) = - mm(3,1)
      nn(1,2) = - xab*mm(1,2)
      nn(1,3) = - xab**2*mm(1,3)
      nn(2,2) =   6.75_R8*sqrt(mab*xab)/(1+xab)**2.5_R8
      nn(2,3) = - 225.0_R8*tab*xab**2/(16.0_R8*(1.0_R8+xab)**3.5_R8)
      nn(3,2) =   nn(2,3)/(xab*tab)
      nn(3,3) =   2625.0_R8*sqrt(tab)*xab**2/
     .           (64.0_R8*(1.0_R8+xab)**4.5_R8)

      return 
      end
