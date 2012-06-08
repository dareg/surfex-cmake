C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FADIES_MT (FA,  KREP, KNUMER, KDATEF )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme permettant d'obtenir la date d'un fichier ouvert
C     pour le logiciel de Fichiers ARPEGE, et deja muni d'une date.
C     ( "DIES" = jour en latin... )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C     (Tableau)  KDATEF (Sortie) ==> Date elle-meme (FA%JPLDAT mots).
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER
      INTEGER KDATEF (FA%JPLDAT)
C
      INTEGER IRANG, J, IREP, INIMES
C
      LOGICAL LLVERF
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FADIES_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ENDIF
C
C         Verrouillage eventuel du fichier.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'ON')
      LLVERF=FA%LFAMUL
C**
C     2.  -  CONTROLE DE DEFINITION PREALABLE DE LA DATE.
C-----------------------------------------------------------------------
C
      IF (FA%LCREAF(IRANG)) THEN
        IREP=-85
        GOTO 1001
      ENDIF
C**
C     3.  -  TRANSFERT DE LA TABLE "FA%MADATE" DANS LE TABLEAU ARGUMENT.
C-----------------------------------------------------------------------
C
      DO 301 J=1,FA%JPLDAT
      KDATEF(J)=FA%MADATE(J,IRANG)
  301 CONTINUE
C
      IREP=0
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
      IF (LLFATA.OR.IXNVMS(IRANG).EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('FADIES_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FADIES'
C
      IF (INIMES.EQ.2) THEN
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', KDATEF(1:5)='',I5,2(''/'',I2),I3,'':'',I2.2,
     S       '', KDATEF(7:8)='',I6,''-'',I6)') KREP,KNUMER,
     S     (KDATEF(J),J=1,5),(KDATEF(J),J=7,8)
        CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S               CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FADIES_MT',1,ZHOOK_HANDLE)
      END

