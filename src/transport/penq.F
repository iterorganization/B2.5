

      subroutine penq(ta,tb,ma,mb,pp,qq)
c--------------------------------------------------------------------
c     subroutine that calculates the coupling coefficients p^ab and 
c     q^ab 
c
c     input    ta   : the temperature of species a
c              tb   : the temperature of species b
c              ma   : the mass of species a
c              mb   : the mass of species b
c     output   pp   : array(2,2) containing the coupling coef. p
c              qq   : array(2,2) containing the coupling coef. q
c--------------------------------------------------------------------
      
      use b2mod_types

      implicit none

      real (kind=R8) ::
     *        ta, tb, ma, mb, pp(2,2), qq(2,2), xab, tab, mab

      tab = ta / tb
      mab = ma / mb
      xab = mab / tab

      pp(1,1) = - 3.*mab / (1+xab)**1.5
      pp(1,2) = -4.5 * xab / (1+xab)**2.5
      pp(2,1) = pp(1,2)
      pp(2,2) = - 3 * xab*(13./4.+ 2*xab + 2.5*xab**2)/
     +            (1+xab)**3.5

      qq(1,1) = - pp(1,1) / tab
      qq(1,2) = - xab*pp(1,2)
      qq(2,1) = - pp(1,2)
      qq(2,2) = 45. * xab**2 / (4.*(1+xab)**3.5)

      return
      end subroutine penq
