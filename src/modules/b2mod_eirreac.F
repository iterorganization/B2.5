      module b2mod_eirreac
      use b2mod_types
      implicit none

      integer nreac
      parameter (nreac=1)
      real (kind=R8) ::
     1 CREAC(9,0:9,NREAC), EARRH(NREAC),
     R RC1MN(NREAC), RC1MX(NREAC), FP1L(NREAC,3), FP1R(NREAC,3),
     R RC2MN(NREAC), RC2MX(NREAC), FP2L(NREAC,3), FP2R(NREAC,3)
      integer
     I ISWR(NREAC), MODCLF(NREAC),
     I IF1MN(NREAC), IF1MX(NREAC), IF2MN(NREAC), IF2MX(NREAC)

      end module b2mod_eirreac
