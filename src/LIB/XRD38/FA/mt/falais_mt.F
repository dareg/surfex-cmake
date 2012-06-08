C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FALAIS_MT (FA,  KREP, KNUMER, CDNOMA, PDONNE, KLONGD )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme de lecture d'un article de donnees non assimila-
C     bles a un champ horizontal sur un fichier ARPEGE.
C       ( Lecture d'un Article Integre Simplement, non code )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                CDNOMA (Entree) ==> Nom de l'article;
C    ( Tableau ) PDONNE (Entree) ==> Donnees a ecrire;
C                KLONGD (Entree) ==> Nombre de mots a ecrire.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER, KLONGD
C
      INTEGER ILCDNO, IREP, IRANG, ILNOMA, INIMES, ILACTI
C
      REAL (KIND=JPDBLR) PDONNE (KLONGD)
C
      LOGICAL LLVERF, LLRLFI
C
      CHARACTER CDNOMA*(*)
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FALAIS_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
      LLRLFI=.FALSE.
      ILCDNO=LEN (CDNOMA)
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ELSEIF (KLONGD.LE.0) THEN
        IREP=-64
        GOTO 1001
      ELSEIF (ILCDNO.LE.0) THEN
        IREP=-65
        GOTO 1001
      ENDIF
C
C         Verrouillage eventuel du fichier.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'ON')
      LLVERF=FA%LFAMUL
C
      IF (FA%LCREAF(IRANG)) THEN
        IREP=-85
        GOTO 1001
      ELSEIF (CDNOMA.EQ.FA%CPCACH.OR.CDNOMA.EQ.FA%CPCADI.OR.
     S    CDNOMA.EQ.FA%CPCAFS.OR.CDNOMA.EQ.FA%CPCARP.OR.
     S    CDNOMA.EQ.FA%CPDATE) THEN
        IREP=-111
        GOTO 1001
      ENDIF
C**
C     2.  -  LECTURE DE L'ARTICLE DE DONNEES SUR LE FICHIER.
C-----------------------------------------------------------------------
C
      ILNOMA=MIN ( FA%NCPCAD, LEN (CDNOMA) )
      CLNOMA(1:ILNOMA)=CDNOMA(1:ILNOMA)
C
      CALL LFILEC_MT (FA%LFI, IREP,KNUMER,CLNOMA(1:ILNOMA),
     S             PDONNE,KLONGD)
      LLRLFI=IREP.NE.0
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
        IF (LHOOK) CALL DR_HOOK('FALAIS_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FALAIS'
C
      IF (IREP.NE.-65) THEN
        ILACTI=MIN (ILCDNO,FA%NCPCAD)
        CLACTI(1:ILACTI)=CDNOMA(:ILACTI)
      ELSE
        ILACTI=8
        CLACTI(1:ILACTI)=FA%CHAINC(:ILACTI)
      ENDIF
C
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', CDNOMA='''''',A,'''''', KLONGD='',I8)')
     S   KREP,KNUMER,CLACTI(1:ILACTI),KLONGD
      CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLACTI(1:ILACTI),LLRLFI)
C
      IF (LHOOK) CALL DR_HOOK('FALAIS_MT',1,ZHOOK_HANDLE)
      END

