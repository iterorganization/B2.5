
      subroutine rswap(r1,r2)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     * r1,r2,r3
      r3=r1
      r1=r2
      r2=r3
      return
      end
