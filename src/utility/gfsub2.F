*//GFSUB2//
*========================================================================
*         S U B R O U T I N E   G F S U B 2
*========================================================================
*
      SUBROUTINE GFSUB2(KARD,NDIMX,NDIMY,NRDX,NRDY,DUMMY)
      use b2mod_types
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
      IMPLICIT NONE
      INTEGER KARD,NDIMX,NDIMY,NRDX,NRDY,
     & LIM,IX,IY,III
      real (kind=R8) ::
     * DUMMY(0:NDIMX+1,0:NDIMY+1)

#ifdef B25_EIRENE
      call read_title(kard)
#endif
      LIM = (NRDX/5)*5 - 4
      DO 110 IY = 1,NRDY
        DO 100 IX = 1,LIM,5
          READ(KARD,910) (DUMMY(IX-1+III,IY),III = 1,5)
  100   CONTINUE
        IF( (LIM+4).EQ.NRDX ) GOTO 110
        READ(KARD,910) (DUMMY(IX,IY),IX = LIM+5,NRDX)
  110 CONTINUE
      RETURN
!pb 27092016   910 FORMAT(5(E16.8))
  910 FORMAT(5(E16.7E3))
*//END GFSUB2//
      END
