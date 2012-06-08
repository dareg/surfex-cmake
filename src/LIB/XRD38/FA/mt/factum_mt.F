C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACTUM_MT (FA,  CDNOMC )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme servant a supprimer un cadre.
C     ( Cadre a TUer Methodiquement ? )
C**
C        Argument : CDNOMC (Entree) ==> Nom symbolique du cadre.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER ILCDNO, IREP, IRANGC, ILNOMC, INIMES, INUMER, J
C
      LOGICAL LLVERG
C
      CHARACTER CDNOMC*(*)
C
C
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  INITIALISATIONS ET CONTROLES SOMMAIRES.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACTUM_MT',0,ZHOOK_HANDLE)
      IF (FA%FACTUM_LLPREA) THEN
C
C         Initialisation eventuelle des variables globales du logiciel.
C
        CALL FARINE_MT (FA, 2)
        FA%FACTUM_LLPREA=.FALSE.
      ENDIF
C
      LLVERG=.FALSE.
      ILCDNO=LEN (CDNOMC)
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
C             Verrouillage global, si necessaire.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
      LLVERG=FA%LFAMUL
C
C          Controle d'existence du cadre specifie.
C
      CALL FANUCA_MT (FA, CDNOMC,IRANGC,.FALSE.)
C
      IF (IRANGC.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ENDIF
C**
C     2.  -  SUPPRESSION PROPREMENT DITE VIA LE SOUS-PROGRAMME "FACTUI".
C-----------------------------------------------------------------------
C
      CALL FACTUI_MT (FA, IREP,IRANGC)
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
      LLFATA=LLMOER(IREP,0)
C
      IF (.NOT.LLFATA.OR.FA%NIMSGA.NE.2)  THEN 
        IF (LHOOK) CALL DR_HOOK('FACTUM_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      INIMES=2
      CLNSPR='FACTUM'
C
      IF (IREP.EQ.-65.AND.ILCDNO.LE.0) THEN
        ILNOMC=8
        CLACTI(1:ILNOMC)=FA%CHAINC(:ILNOMC)
      ELSE
        ILNOMC=MIN (LEN (CLACTI),ILNOMC)
        CLACTI=CDNOMC(1:ILNOMC)
      ENDIF
C
      WRITE (UNIT=CLMESS,FMT='(''CDNOMC='''''',A,'''''''')')
     S     CLACTI(1:ILNOMC)
      INUMER=FA%JPNIIL
      CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR, CLACTI(1:ILNOMC),.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FACTUM_MT',1,ZHOOK_HANDLE)
      END

