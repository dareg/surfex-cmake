C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACTUI_MT (FA,  KREP, KRANGC )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme A USAGE INTERNE AU LOGICIEL. Fait la suppression
C     d'un cadre ( vis-a-vis des tables du logiciel ) .
C        En mode multi-taches, il doit y avoir verrouillage global
C     de la zone d'appel au sous-programme.
C**
C        Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                    KRANGC (Entree) ==> Rang du cadre dans les tables.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANGC
C
      INTEGER J, IPOSCA, INIMES, INUMER
#include "facom_mt.h"
C**
C     1.  -  CONTROLE PREALABLE DE COHERENCE.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACTUI_MT',0,ZHOOK_HANDLE)
      IF (KRANGC.LE.0.OR.KRANGC.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C**
C     2.  -  RECHERCHE DU CADRE DANS LA TABLE "FA%NCAIND".
C-----------------------------------------------------------------------
C
      DO 201 J=1,FA%NCADEF
C
      IF (FA%NCAIND(J).EQ.KRANGC) THEN
        IPOSCA=J
        GOTO 202
      ENDIF
C
  201 CONTINUE
C
      KREP=-66
      GOTO 1001
C
  202 CONTINUE
C
      IF (FA%NULCAD(IPOSCA).NE.0) THEN
        KREP=-67
        GOTO 1001
      ENDIF
C**
C     3.  -  MISE A JOUR DES TABLES.
C-----------------------------------------------------------------------
C
      FA%CNOMCA(KRANGC)=' '
      FA%NCADEF=FA%NCADEF-1
C
      DO 301 J=IPOSCA,FA%NCADEF
      FA%NCAIND(J)=FA%NCAIND(J+1)
  301 CONTINUE
C
      KREP=0
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE EVENTUELLE,
C            VIA LE sous-programme "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
      LLFATA=KREP.EQ.-66.OR.(KREP.NE.0.AND.FA%NRFAGA.NE.2)
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FACTUI'
        INUMER=FA%JPNIIL
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANGC='',I4)')
     S     KREP,KRANGC
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S               CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FACTUI_MT',1,ZHOOK_HANDLE)
      END

