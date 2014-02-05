
*//GFSUB3R//
*=======================================================================
*                 S U B R O U T I N E  G F S U B 3
*=======================================================================
*
*
*  1. PURPOSE
*
*     GFSUB3R WRITES AN ARRAY ONTO AN OUTPUT FILE
*
*
*=======================================================================
*//DECLARATIONS//

        SUBROUTINE GFSUB3R(KARD,NX,NY,NDIMX,NDIMY,NDIMS,DUMMY)
      use b2mod_types
        IMPLICIT NONE
        INTEGER KARD,NX,NY,NDIMX,NDIMY, NDIMS, ND1, LIM, IS, IY, IX, III
        REAL (KIND=R8) ::
     *   DUMMY(0:NDIMX+1,0:NDIMY+1,1:NDIMS)
        ND1 = NX + 2
        LIM = (ND1/5)*5 - 4
        do  120  IS = 1,NDIMS
        do  110  IY = -1,NY
        do  100  IX = 1,LIM,5
        READ(KARD,910) (DUMMY(IX-2+III,IY,IS),III = 1,5)
  100   CONTINUE
        IF( (LIM+4) .EQ. ND1 )   GOTO 110
        READ(KARD,910) (DUMMY(IX-1,IY,IS),IX = LIM+5,ND1)
  110   CONTINUE
  120   CONTINUE
        RETURN
  910   FORMAT(5(E16.8))
*=======================================================================
*//END GFSUB3R//

        END
