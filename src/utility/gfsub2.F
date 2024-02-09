*//GFSUB2//
*========================================================================
*         S U B R O U T I N E   G F S U B 2
*========================================================================
*
      SUBROUTINE GFSUB2(KARD,NDIMX,NRDX,DUMMY)
      use b2mod_types
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
      IMPLICIT NONE
      INTEGER KARD,NDIMX,NRDX,
     & LIM,IX,III
      real (kind=R8) ::
     * DUMMY(NDIMX)

#ifdef B25_EIRENE
      call read_title(kard)
#endif
      LIM = (NRDX/5)*5 - 4
      DO 100 IX = 1,LIM,5
        READ(KARD,910) (DUMMY(IX-1+III),III = 1,5)
  100 CONTINUE
      IF( (LIM+4).NE.NRDX )
     . READ(KARD,910) (DUMMY(IX),IX = LIM+5,NRDX)
      RETURN
!pb 27092016   910 FORMAT(5(E16.8))
  910 FORMAT(5(E16.7E3))
*//END GFSUB2//
      END SUBROUTINE GFSUB2
