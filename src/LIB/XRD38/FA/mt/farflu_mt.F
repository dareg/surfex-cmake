C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FARFLU_MT (FA,  KXNIVV, KXTRON, KXLATI, KXLONG )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme servant a specifier des Limites Utilisateur
C     en termes de Resolutions horizontale et verticale, valides
C     globalement pour tous les Fichiers et Cadres ARPEGE.
C            ( Resolution Fichiers - Limites Utilisateur )
C**
C        Arguments : KXNIVV ==> Nombre maximum de niveaux verticaux;
C  (tous d'Entree)   KXTRON ==> Troncature maximum;
C                    KXLATI ==> Nombre maximum de latitudes pole a pole;
C                    KXLONG ==> Nombre maxi de longitudes par parallele.
C*
C        S'il y a des cadres deja definis (dynamiquement ou non) avec
C     des valeurs correspondantes plus grandes, cela provoque une erreur
C     globale.
C        Si les valeurs donnees en argument depassent les limites cor-
C     repondantes du logiciel, elles sont ecretees, avec un message.
C        Une messagerie de niveau 1 est emise dans le cas normal ou
C     ci-dessus.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KXNIVV, KXTRON, KXLATI, KXLONG
C
      INTEGER INUMER, INIMES, IREP, IMINIM, IXNIVV, IXTRON, IXLATI, J
      INTEGER IXLONG, IRANGC
C
      LOGICAL LLDEPA, LLVERG
C
C
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  SI PREMIERE UTILISATION, APPEL AU SOUS-PROGRAMME "FARINE".
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FARFLU_MT',0,ZHOOK_HANDLE)
      IF (FA%FARFLU_LLPREA) THEN
        CALL FARINE_MT (FA, 2)
        FA%FARFLU_LLPREA=.FALSE.
      ENDIF
C**
C     2.  -  CONTROLES.
C-----------------------------------------------------------------------
C
      IMINIM=MIN (KXNIVV,KXTRON,KXLATI,KXLONG)
C
      IF (IMINIM.LE.0) THEN
        LLVERG=.FALSE.
        LLDEPA=.FALSE.
        IREP=-64
        GOTO 1001
      ELSE
        LLDEPA=KXNIVV.GT.FA%JPXNIV.OR.KXTRON.GT.FA%JPXTRO.OR.
     S         KXLATI.GT.FA%JPXLAT.OR.KXLONG.GT.FA%JPXLON
        IXNIVV=MIN (KXNIVV,FA%JPXNIV)
        IXTRON=MIN (KXTRON,FA%JPXTRO)
        IXLATI=MIN (KXLATI,FA%JPXLAT)
        IXLONG=MIN (KXLONG,FA%JPXLON)
      ENDIF
C
C             Verrouillage global, si necessaire.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
      LLVERG=FA%LFAMUL
C
      IF (IXTRON.LE.FA%NSTROI) THEN
        IREP=-99
        GOTO 1001
      ENDIF
C
C        Controles vis-a-vis des cadres deja definis, s'il y en a.
C
      DO 201 J=1,FA%NCADEF
      IRANGC=FA%NCAIND(J)
C
      IF (FA%MTRONC(IRANGC).GT.IXTRON) THEN
        IREP=-104
        GOTO 1001
      ELSEIF (FA%NNIVER(IRANGC).GT.IXNIVV) THEN
        IREP=-105
        GOTO 1001
      ELSEIF (FA%NLATIT(IRANGC).GT.IXLATI) THEN
        IREP=-106
        GOTO 1001
      ELSEIF (FA%NXLOPA(IRANGC).GT.IXLONG) THEN
        IREP=-107
        GOTO 1001
      ENDIF
C
  201 CONTINUE
C**
C     3.  -  MODIFICATION DES LIMITES USAGER.
C-----------------------------------------------------------------------
C
      FA%NXNIVV=IXNIVV
      FA%NXTRON=IXTRON
      FA%NXLATI=IXLATI
      FA%NXLONG=IXLONG
      IREP=0
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE sous-programme "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
C          Deverrouillage global eventuel.
C
      IF (LLVERG) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'OFF')
C
      LLFATA=IREP.NE.0.AND.FA%NRFAGA.NE.2
C
      IF (LLFATA) THEN
        INIMES=2
      ELSEIF (FA%NIMSGA.EQ.0) THEN
        IF (LHOOK) CALL DR_HOOK('FARFLU_MT',1,ZHOOK_HANDLE)
        RETURN
      ELSE
        INIMES=FA%NIMSGA
      ENDIF
C
      CLNSPR='FARFLU'
      INUMER=FA%JPNIIL
C
      IF (INIMES.EQ.2) THEN
        WRITE (UNIT=CLMESS,FMT='(''KXNIVV='',I4,'', KXTRON='',I4,
     S         '', KXLATI='',I4,'', KXLONG='',I4)')
     S     KXNIVV,KXTRON,KXLATI,KXLONG
        CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S               CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (INIMES.GE.1) THEN
C
        IF (LLDEPA) THEN
          WRITE (UNIT=CLMESS,FMT=
     S   '(''LIMITES USAGER (***ECRETEES***): FA%NXNIVV='',I4,
     S   '', FA%NXTRON='',I4,'', FA%NXLATI='',I4,'', FA%NXLONG='',I4)')
     S     FA%NXNIVV,FA%NXTRON,FA%NXLATI,FA%NXLONG
        ELSE
          WRITE (UNIT=CLMESS,FMT=
     S    '(''LIMITES USAGER EFFECTIVES: FA%NXNIVV='',I4,
     S    '', FA%NXTRON='',I4,'', FA%NXLATI='',I4,'', FA%NXLONG='',I4)')
     S     FA%NXNIVV,FA%NXTRON,FA%NXLATI,FA%NXLONG
        ENDIF
C
        CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,.FALSE.,CLMESS,
     S               CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FARFLU_MT',1,ZHOOK_HANDLE)
      END

