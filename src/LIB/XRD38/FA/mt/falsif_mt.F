C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FALSIF_MT (FA,  KREP, KNUMER, CDIDEN )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme renvoyant le NOM de l'Identificateur
C     d'un fichier ARPEGE.
C       ( Lecture Specifique de l'Identificateur de Fichier )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                CDIDEN (Sortie) ==> Nom de l'identificateur.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER
C
      INTEGER IREP, ILIDEN, IRANG, J, ILONGN, INIMES, ILACTI
C
      LOGICAL LLVERF, LLRLFI
C
      CHARACTER CDIDEN*(*)
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FALSIF_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
      LLRLFI=.FALSE.
      ILIDEN=LEN (CDIDEN)
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ELSEIF (ILIDEN.LE.0) THEN
        IREP=-65
        GOTO 1001
      ELSE
        IREP=0
      ENDIF
C
C         Verrouillage eventuel du fichier.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'ON')
      LLVERF=FA%LFAMUL
C**
C     2.  -  ON RENVOIE LE NOM D'IDENTIFICATEUR, APRES CONTROLE EVENTUEL
C            D'UNE VARIABLE CARACTERE SUFISAMMENT LONGUE.
C-----------------------------------------------------------------------
C
      IF (ILIDEN.GE.FA%NCPCAD) THEN
        CDIDEN=FA%CIDENT(IRANG)
      ELSE
C
        DO 201 J=FA%NCPCAD,1,-1
C
        IF (FA%CIDENT(IRANG)(J:J).NE.' ') THEN
          ILONGN=J
          GOTO 202
        ENDIF
C
  201   CONTINUE
C
        IREP=-66
        GOTO 1001
C
  202   CONTINUE
C
        IF (ILONGN.GT.ILIDEN) THEN
          IREP=-69
        ELSE
          CDIDEN=FA%CIDENT(IRANG)(1:ILONGN)
        ENDIF
C
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG)
C
C        Deverrouillage eventuel du fichier.
C
      IF (LLVERF) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'OFF')
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=IXNVMS(IRANG)
      ENDIF
C
      IF (.NOT.LLFATA.AND.INIMES.NE.2)  THEN 
        IF (LHOOK) CALL DR_HOOK('FALSIF_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FALSIF'
C
      IF (IREP.EQ.0.OR.IREP.EQ.-69) THEN
        ILACTI=FA%NCPCAD
        CLACTI=FA%CIDENT(IRANG)(1:ILACTI)
      ELSE
        ILACTI=8
        CLACTI(1:ILACTI)=FA%CHAINC(:ILACTI)
      ENDIF
C
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', CDIDEN='''''',A,'''''''')')
     S   KREP,KNUMER,CLACTI(1:ILACTI)
      CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLACTI(1:ILACTI),LLRLFI)
C
      IF (LHOOK) CALL DR_HOOK('FALSIF_MT',1,ZHOOK_HANDLE)
      RETURN
      END

