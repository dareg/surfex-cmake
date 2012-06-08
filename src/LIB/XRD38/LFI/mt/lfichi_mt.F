C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFICHI_MT (LFI, KREP, CDSTRU, KVAL, KPOSC2 )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        SOUS-PROGRAMME *INTERNE* DU LOGICIEL DE FICHIERS INDEXES LFI;
C     permet de decoder une valeur entiere (CHIffres)
C     dans une chaine de caracteres.
C**
C    ARGUMENTS : KREP   (Sortie) ==> Code-Reponse du sous-programme;
C                CDSTRU (Entree) ==> Chaine a decoder;
C                KVAL   (Sortie) ==> Valeur entiere decodee;
C                KPOSC2 (Sortie) ==> Position du dernier chiffre.
C
      TYPE(LFICOM) :: LFI
      CHARACTER CDSTRU*(*), CLFORM*7
C
      INTEGER KREP, KVAL, KPOSC2
      INTEGER ILUSTR, J, IPOSC1, IPOSC2
#include "lficom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFICHI_MT',0,ZHOOK_HANDLE)
      KREP=0
      ILUSTR=LEN (CDSTRU)
C
      DO 222 J=1,ILUSTR
C
      IF (CDSTRU(J:J).NE.' ') THEN
C
        IPOSC1=INDEX (LFI%LFICHI_CLCHIF,CDSTRU(J:J))
C
        IF (IPOSC1.EQ.0) THEN
          KREP=-40
          GOTO 1001
        ENDIF
C
        IPOSC1=J
        GOTO 223
C
      ENDIF
C
  222 CONTINUE
C
      IPOSC1=1
C
  223 CONTINUE
C
      DO 224 J=IPOSC1+1,ILUSTR
C
      IF (INDEX (LFI%LFICHI_CLCHIF,CDSTRU(J:J)).EQ.0) THEN
        IPOSC2=J-1
        GOTO 225
      ENDIF
C
  224 CONTINUE
C
      IPOSC2=ILUSTR
C
  225 CONTINUE
C
      WRITE (UNIT=CLFORM,FMT='(''(BN,I'',I1,'')'')') IPOSC2-IPOSC1+1
      READ (UNIT=CDSTRU(IPOSC1:IPOSC2),FMT=CLFORM,ERR=226) KVAL
      KPOSC2=IPOSC2
      GOTO 1001
C
  226 CONTINUE
C
      KREP=-40
C
 1001 CONTINUE
C
      IF (LHOOK) CALL DR_HOOK('LFICHI_MT',1,ZHOOK_HANDLE)
      END

