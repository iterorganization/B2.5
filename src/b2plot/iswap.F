
      subroutine iswap(i1,i2)
      use b2mod_types
      Implicit None
      Integer :: i1, i2, i3
      i3=i1
      i1=i2
      i2=i3
      return
      end
