      integer nreac
      parameter (nreac=1)
      real (kind=R8) ::
     1 CREAC(9,0:9,NREAC),
     R RCMN(NREAC),RCMX(NREAC),FPARM(NREAC,6)
      integer
     I ISWR(NREAC),MODCLF(NREAC),IFEXMN(NREAC),IFEXMX(NREAC)
      COMMON /COMAMF_b2plot/
C  NAMF, REAL
     R CREAC,RCMN,RCMX,FPARM,
C  MAMF, INTEGER
     I ISWR,MODCLF,IFEXMN,IFEXMX
