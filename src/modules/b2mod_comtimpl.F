      module b2mod_comtimpl
      use b2mod_types
      implicit none
 
#ifndef USE_CDF_B2STIME
C-----------------------------------------------------------------------
C
CCTIMPL      COMMON FOR PLOTS OF VARIABLES OVER TIME
C
      INTEGER    NTDIM
      PARAMETER (NTDIM = 10000)
C
C     W A R N U N G: WENN NTDIM GEANDERT WIRD, BITTE AUCH DIE
C     *************  DIMENSIONEN DER HILFSARRAYS IN "B2DRIV" AENDERN!
C
      INTEGER  ITSTEP, NTSTEP, ILAST, NKK1, NKK2
      REAL (kind=R8) ::
     *  NE3DMX, T3DLMX, T3DIMX, T3DRMX, TIMES
      REAL (kind=R8), allocatable ::
     *  NE3DL(:,:), NE3DI(:,:), NE3DR(:,:), TE3DL(:,:), TE3DI(:,:), 
     *  TE3DR(:,:),
     *  TI3DL(:,:), TI3DI(:,:), TI3DR(:,:), NE3DA(:,:), TE3DA(:,:),  
     *  TI3DA(:,:), AJ(:,:),
     *  FNIXIP(:), FEEXIP(:), FEIXIP(:), FNIXAP(:), FEEXAP(:),  
     *  FEIXAP(:),
     *  FNISIP(:), FEESIP(:), FEISIP(:), FNISAP(:), FEESAP(:),  
     *  FEISAP(:),
     *  FNISIPP(:), FEESIPP(:), FEISIPP(:), FNISAPP(:),  
     *  FEESAPP(:), FEISAPP(:),
     *  FNIYIP(:), FEEYIP(:), FEIYIP(:), FNIYAP(:), FEEYAP(:),  
     *  FEIYAP(:),
     *  NESEPI(:), TESEPI(:), TISEPI(:), NESEPM(:), TESEPM(:),  
     *  TISEPM(:),
     *  NESEPA(:), TESEPA(:), TISEPA(:),
     *  NEMXIP(:), TEMXIP(:), TIMXIP(:), NEMXAP(:), TEMXAP(:),  
     *  TIMXAP(:), TIMESA(:),
     *  fn3dl(:,:),fn3dr(:,:),fe3dl(:,:),fe3dr(:,:),fi3dl(:,:), 
     *  fi3dr(:,:),PWMXIP(:),PWMXAP(:),
     1  tmne(:),tmte(:),tmti(:),tmhacore(:),tmhasol(:),tmhadiv(:),
     2  an3dl(:,:),an3di(:,:),an3da(:,:),an3dr(:,:),mn3dl(:,:), 
     *  mn3di(:,:),mn3da(:,:),mn3dr(:,:)
      LOGICAL  LSTORE

      contains
      
      subroutine allocate_comtimpl(ny)
      integer ny

      allocate(
     >      NE3DL(0:NY+1,NTDIM)   , NE3DI(0:NY+1,NTDIM)  ,
     >      NE3DA(0:NY+1,NTDIM)   , NE3DR(0:NY+1,NTDIM)  ,
     >      TE3DL(0:NY+1,NTDIM)   , TI3DL(0:NY+1,NTDIM)  ,
     >      TE3DI(0:NY+1,NTDIM)   , TI3DI(0:NY+1,NTDIM)  ,
     >      TE3DA(0:NY+1,NTDIM)   , TI3DA(0:NY+1,NTDIM)  ,
     >      TE3DR(0:NY+1,NTDIM)   , TI3DR(0:NY+1,NTDIM)  ,
     >      TIMESA(NTDIM)         ,    AJ(0:NY+1,NTDIM)  ,
     >      fn3dl(0:NY+1,NTDIM)   , fn3dr(0:NY+1,NTDIM)  ,
     >      fe3dl(0:NY+1,NTDIM)   , fe3dr(0:NY+1,NTDIM)  ,
     >      fi3dl(0:NY+1,NTDIM)   , fi3dr(0:NY+1,NTDIM)  ,
     >      an3dl(0:NY+1,NTDIM)   , an3di(0:NY+1,NTDIM)  ,
     >      an3da(0:NY+1,NTDIM)   , an3dr(0:NY+1,NTDIM)  ,
     >      mn3dl(0:NY+1,NTDIM)   , mn3di(0:NY+1,NTDIM)  ,
     >      mn3da(0:NY+1,NTDIM)   , mn3dr(0:NY+1,NTDIM)  ,
     >      FNIXIP(NTDIM)   , FEEXIP(NTDIM)   , FEIXIP(NTDIM)   ,
     >      FNIXAP(NTDIM)   , FEEXAP(NTDIM)   , FEIXAP(NTDIM)   ,
     >      FNISIP(NTDIM)   , FEESIP(NTDIM)   , FEISIP(NTDIM)   ,
     >      FNISAP(NTDIM)   , FEESAP(NTDIM)   , FEISAP(NTDIM)   ,
     >      FNISIPP(NTDIM)  , FEESIPP(NTDIM)  , FEISIPP(NTDIM)  ,
     >	    FNISAPP(NTDIM)  , FEESAPP(NTDIM)  , FEISAPP(NTDIM)  ,
     >      FNIYIP(NTDIM)   , FEEYIP(NTDIM)   , FEIYIP(NTDIM)   ,
     >      FNIYAP(NTDIM)   , FEEYAP(NTDIM)   , FEIYAP(NTDIM)   ,
     >      NESEPI(NTDIM)   , TESEPI(NTDIM)   , TISEPI(NTDIM)   ,
     >      NESEPM(NTDIM)   , TESEPM(NTDIM)   , TISEPM(NTDIM)   ,
     >      NESEPA(NTDIM)   , TESEPA(NTDIM)   , TISEPA(NTDIM)   ,
     >      NEMXIP(NTDIM)   , TEMXIP(NTDIM)   , TIMXIP(NTDIM)   ,
     >      NEMXAP(NTDIM)   , TEMXAP(NTDIM)   , TIMXAP(NTDIM)   ,
     >      PWMXIP(NTDIM)   , PWMXAP(NTDIM)   ,
     >	    tmne(ntdim),tmte(ntdim),tmti(ntdim),
     >      tmhacore(ntdim),tmhasol(ntdim),tmhadiv(ntdim))

      end subroutine allocate_comtimpl
C
#else
C-----------------------------------------------------------------------
C
CCTIMPL      COMMON FOR PLOTS OF VARIABLES OVER TIME
C
      INTEGER    NTDIM
      PARAMETER (NTDIM = 100)
      INTEGER  NTSTEP
C
#endif
      end module b2mod_comtimpl
