C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAUTIF_MT (FA,  KREP, KNUMER, CDIDEN )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme permettant de donner un NOM a l'Identificateur
C     d'un fichier ARPEGE.
C       ( l'Utilisateur Traite son Identificateur de Fichier )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                CDIDEN (Entree) ==> Nom de l'identificateur.
C
C       Une messagerie de niveau 1 est emise dans les cas "normaux"
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER
C
      INTEGER IREP, ILIDEN, IRANG, INIMES, ILACTI
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
      IF (LHOOK) CALL DR_HOOK('FAUTIF_MT',0,ZHOOK_HANDLE)
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
      ENDIF
C
C         Verrouillage eventuel du fichier.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'ON')
      LLVERF=FA%LFAMUL
      CLNOMA=FA%CIDENT(IRANG)
C
      IF (CDIDEN.EQ.FA%CPCACH.OR.CDIDEN.EQ.FA%CPCADI.OR.
     S    CDIDEN.EQ.FA%CPCAFS.OR.CDIDEN.EQ.FA%CPCARP.OR.
     S    CDIDEN.EQ.FA%CPDATE) THEN
        IREP=-111
        GOTO 1001
      ENDIF
C**
C     2.  -  ON RENOMME L'ARTICLE IDENTIFICATEUR, QUI EXISTE TOUJOURS SI
C            LE FICHIER EST OUVERT, AU MOINS AVEC UN NOM PAR DEFAUT.
C-----------------------------------------------------------------------
C
      IF (CDIDEN.NE.FA%CIDENT(IRANG)) THEN
        CALL LFIREN_MT (FA%LFI, IREP,KNUMER,FA%CIDENT(IRANG),CDIDEN)
        LLRLFI=IREP.NE.0
        IF (.NOT.LLRLFI) FA%CIDENT(IRANG)=CDIDEN
      ELSE
        IREP=0
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
      IF (.NOT.LLFATA.AND.INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('FAUTIF_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FAUTIF'
C
      IF (IREP.NE.-65) THEN
        ILACTI=FA%NCPCAD
        CLACTI(1:ILACTI)=CDIDEN(1:MIN (ILIDEN,ILACTI))
      ELSE
        ILACTI=8
        CLACTI(1:ILACTI)=FA%CHAINC(:ILACTI)
      ENDIF
C
      IF (INIMES.EQ.2) THEN
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S         '', CDIDEN='''''',A,'''''''')')
     S   KREP,KNUMER,CLACTI(1:ILACTI)
        CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,
     S               CLMESS,CLNSPR,
     S               CLACTI(1:ILACTI),LLRLFI)
      ENDIF
C
C        La messagerie qui suit n'est pas emise en cas d'erreur fatale.
C
      IF (INIMES.GE.1.AND.IRANG.NE.0) THEN
        WRITE (UNIT=CLMESS,FMT=
     S  '(''Ancien Identificateur de l''''unite logique'',I3,
     S    '' : '''''',A,'''''', Nouveau: '''''',A,'''''''')')
     S  KNUMER,CLNOMA,FA%CIDENT(IRANG)
        CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI(1:ILACTI),.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAUTIF_MT',1,ZHOOK_HANDLE)
      END

