
      SUBROUTINE COLOR2(N)
      use b2mod_types
      integer n
      integer i, iw, icnt
      real hues, xhue, red, green, blue
#include <graphics_vars.h>
      do iw=1,nwork
C BACKGROUND COLOR
C WHITE
          CALL GSCR(iw,0,.6*1.,.6*1.,.6*1.)
C BLACK
          CALL GSCR(iw,1,0.,0.,0.)
        enddo
C Choose other foreground colors spaced equally around the spectrum
        ICNT=0
        HUES=360./(N+2)
        DO I=1,N
          XHUE=(n-i)*HUES
          CALL HLSRGB(XHUE,60.,75.,RED,GREEN,BLUE)
          IF (XHUE.LE.36.0) THEN
            do iw=1,nwork
              CALL GSCR(iw,N+3-I,RED,GREEN,BLUE)
            enddo
            ICNT=ICNT+1
          ELSE
              do iw=1,nwork
                CALL GSCR(iw,I-ICNT+2,RED,GREEN,BLUE)
              enddo
            ENDIF
          ENDDO
          RETURN
        END
