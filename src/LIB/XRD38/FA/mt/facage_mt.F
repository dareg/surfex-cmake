C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACAGE_MT (FA,  CDNOMC, LDGARD )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme servant a redefinir l'option de conservation
C     d'un cadre preexistant ( CAdre a Garder Eventuellement... )
C**
C        Arguments : CDNOMC ==> Nom symbolique du cadre;
C  (tous d'Entree)   LDGARD ==> Vrai si le cadre doit etre conserve meme
C                               apres la fermeture du dernier fichier
C                               qui s'y rattache.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER ILCDNO, ILNOMC, J, IRANGC, IREP, INIMES, INUMER
C
      LOGICAL LLVERG, LDGARD
C
      CHARACTER CDNOMC*(*)
#include "facom_mt.h"
C
C
C**
C     0.  -  SI PREMIERE UTILISATION, APPEL AU SOUS-PROGRAMME "FARINE".
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACAGE_MT',0,ZHOOK_HANDLE)
      IF (FA%FACAGE_LLPREA) THEN
        CALL FARINE_MT (FA, 2)
        FA%FACAGE_LLPREA=.FALSE.
      ENDIF
C**
C     1.  -  CONTROLE DE L'ARGUMENT "CDNOMC".
C-----------------------------------------------------------------------
C
      LLVERG=.FALSE.
      ILCDNO=LEN (CDNOMC)
      ILNOMC=1
C
      IF (ILCDNO.LE.0) THEN
        IREP=-65
        GOTO 1001
      ELSEIF (CDNOMC.EQ.' ') THEN
        IREP=-68
        GOTO 1001
      ENDIF
C
      DO 101 J=ILCDNO,1,-1
C
      IF (CDNOMC(J:J).NE.' ') THEN
        ILNOMC=J
        GOTO 102
      ENDIF
C
  101 CONTINUE
C
  102 CONTINUE
C
      IF (ILNOMC.GT.FA%NCPCAD) THEN
        IREP=-65
        GOTO 1001
      ENDIF
C**
C     2.  -  RECHERCHE DU CADRE DANS LES TABLES.
C-----------------------------------------------------------------------
C
C             Verrouillage global prealable, si necessaire.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
      LLVERG=FA%LFAMUL
C
      CALL FANUCA_MT (FA, CDNOMC,IRANGC,.FALSE.)
C
      IF (IRANGC.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ENDIF
C**
C     3.  -  MISE A JOUR DU NIVEAU DE CONSERVATION.
C-----------------------------------------------------------------------
C
      IF (LDGARD) THEN
        FA%NGARDE(IRANGC)=2
      ELSE
        FA%NGARDE(IRANGC)=0
      ENDIF
C
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
      IF (LLFATA.OR.FA%NIMSGA.EQ.2) THEN
        INIMES=2
        CLNSPR='FACAGE'
C
        IF (IREP.EQ.-65.AND.ILCDNO.LE.0) THEN
          ILNOMC=8
          CLACTI(1:ILNOMC)=FA%CHAINC(:ILNOMC)
        ELSE
          ILNOMC=MIN (LEN (CLACTI),ILNOMC)
          CLACTI(1:ILNOMC)=CDNOMC(1:ILNOMC)
        ENDIF
C
        ILNOMC=MIN (ILNOMC,FA%NCPCAD)
        WRITE (UNIT=CLMESS,
     S         FMT='(''CDNOMC= '''''',A,'''''', LDGARD= '',
     S         L1,'', CODE INTERNE='',I4)') 
     S         CLACTI(1:ILNOMC),LDGARD,IREP
        INUMER=FA%JPNIIL
        CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S               CLNSPR,CLACTI(1:ILNOMC),.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FACAGE_MT',1,ZHOOK_HANDLE)
      END

