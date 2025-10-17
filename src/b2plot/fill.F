
C ================= User-defined routine for contour filling ======
      SUBROUTINE FILL(XWRK,YWRK,NWRK,IAREA,IGRP,NGRPS)
C =================================================================
      use b2mod_types
      use b2mod_graphics
      implicit none
      INTEGER NWRK, IAREA, IGRP, NGRPS
      REAL XWRK, YWRK
      DIMENSION XWRK(*),YWRK(*),IAREA(*),IGRP(*)
      LOGICAL DEBUG,PAINT
      INTEGER IDFILL, IGR1, IGR2, IAR1, IAR2, ISWP
      integer choose_color

#ifdef DBG
      DEBUG=.TRUE.
#else
      DEBUG=.FALSE.
#endif
      IDFILL=-1
      PAINT=.FALSE.
      IGR1=IGRP(1)
      IAR1=IAREA(1)
      IGR2=IGRP(2)
      IAR2=IAREA(2)
      IF(IGR1.NE.1) THEN
        ISWP=IGR1
        IGR1=IGR2
        IGR2=ISWP
        ISWP=IAR1
        IAR1=IAR2
        IAR2=ISWP
C       WRITE (*,*) 'swapped 1 and 2 ...'
      ENDIF
      IF(IAR1.LT.1) THEN
        RETURN
      ELSEIF(IAR2.LT.1) THEN
        RETURN
      ELSE
        IF(IAR2.GE.1) THEN
          PAINT=.TRUE.
          IDFILL=IAR1
        ELSE
          RETURN
        ENDIF
      ENDIF
      IF(PAINT) THEN
        IF(DEBUG) WRITE (*,*) 'painted IGRP',IGR1,'in IAREA',IDFILL
        CALL GSFACI(choose_color(IDFILL-1,0))
        CALL GFA(NWRK-1,XWRK,YWRK)
      ENDIF
C Otherwise, do nothing
      RETURN
      END SUBROUTINE FILL
