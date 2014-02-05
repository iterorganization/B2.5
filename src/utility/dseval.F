      FUNCTION DSEVAL(N,U,X,Y,B,C,D)
      use b2mod_types
      implicit none
      INTEGER N
      REAL (KIND=R8) :: U,X(N),Y(N),B(N),C(N),D(N), DSEVAL
C
C  ÝTA ÏOÄÏPOÃPAMMA BÛ×ÈCËßET ÇHA×EHÈE KÓÁÈ×ECKOÃO
C  CÏËAÉHA
C
C  DSEVAL=Y(I)+B(I)*(U-X(I))+C(I)*(U-X(I))**2+D(I)*(U-X(I))**3
C
C  ÃÄE X(I).LT.U.LT.X(I+1). ÈCÏOËÜÇÓETCß CXEMA ÃOPHEPA
C
C  ECËÈ U.LT.X(1), TO ÁEPETCß ÇHA×EHÈE I=1.
C  ECËÈ U.GE.X(N), TO ÁEPETCß ÇHA×EHÈE I=N.
C
C  BXOÄHAß ÈHÔOPMAÖÈß.
C
C     N     -×ÈCËO ÇAÄAHHÛX TO×EK
C     U     -AÁCÖÈCCA, ÄËß KOTOPOÉ BÛ×ÈCËßETCß ÇHA×EHÈE
C            CÏËAÉHA
C     X,Y   -MACCÈBÛ ÇAÄAHHÛX AÁCÖÈCC È OPÄÈHAT
C     B,C,D -MACCÈBÛ KOÝÔÔÈÖÈEHTOB CÏËAÉHA, BÛ×ÈCËEHHÛE
C            ÏOÄÏPOÃPAMMOÉ SPLINE
C
C  ECËÈ ÏO CPABHEHÈÞ C ÏPEÄÛÄÓÙÈM BÛÇOBOM U HE
C  HAXOÄÈTCß B TOM ÆE ÈHTEPBAËE, TO ÄËß PAÇÛCKAHÈß
C  HÓÆHOÃO ÈHTEPBAËA ÏPÈMEHßETCß ÄBOÈ×HÛÉ ÏOÈCK.
C
      INTEGER I,J,K
      REAL (KIND=R8) :: DX
	save i
      DATA I/1/
      IF(I.GE.N) I=1
      IF(U.LT.X(I)) GO TO 10
      IF(U.LE.X(I+1)) GO TO 30
C
C  ÄBOÈ×HÛÉ ÏOÈCK
C
 10   I=1
      J=N+1
 20   K=(I+J)/2
      IF(U.LT.X(K))J=K
      IF(U.GE.X(K))I=K
      IF(J.GT.I+1)GO TO 20
C
C  BÛ×ÈCËÈTÜ CÏËAÉH
C
 30   DX=U-X(I)
      DSEVAL=Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
      RETURN
      END
