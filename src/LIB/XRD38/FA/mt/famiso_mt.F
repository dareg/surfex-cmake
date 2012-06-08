C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAMISO_MT (FA,  LDEBUG )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'activer ou de desactiver le mode
C     "Mise au point du logiciel". ( par defaut, inactif )
C        A noter que le mode "mise au point" du logiciel LFI n'est pas
C     modifie.
C**
C        Argument : LDEBUG (Entree) ==> Vrai si on doit activer ce mode.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER IREP, INUMER, INIMES
C
      LOGICAL LDEBUG
#include "facom_mt.h"
C
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAMISO_MT',0,ZHOOK_HANDLE)
      IF (FA%FAMISO_LLPREA) THEN
        CALL FARINE_MT (FA, 2)
        FA%FAMISO_LLPREA=.FALSE.
      ENDIF
C
      FA%LFAMOP=LDEBUG
C
C  Prise en compte du niveau de messagerie dans GRIBEX
C
      IF (FA%LFAMOP) THEN
        CALL GRSDBG(1)
      ELSE
        CALL GRSDBG(0)
      ENDIF
C
C        MESSAGERIE EVENTUELLE .
C
      IF (FA%NIMSGA.NE.2)  THEN 
        IF (LHOOK) CALL DR_HOOK('FAMISO_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      INUMER=FA%JPNIIL
      INIMES=2
      IREP=0
      CLNSPR='FAMISO'
      WRITE (UNIT=CLMESS,FMT='(''LDEBUG= '',L1)') LDEBUG
      CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,.FALSE.,CLMESS,
     S                CLNSPR,CLACTI,.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FAMISO_MT',1,ZHOOK_HANDLE)
      END

