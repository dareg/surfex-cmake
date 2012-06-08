C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIFRA_MT (LFI, LDFRAN )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce Sous-Programme permet aux messages (ulterieurs) emis par le
C     logiciel LFI d'etre rediges en FRANCAIS ou en ANGLAIS.
C        L'option par defaut est definie dans le sous-programme *LFIINI*
C
C        This Subroutine enables (further) messages from LFI software
C     to be written in FRENCH or ENGLISH. Default mode is defined in
C     routine *LFIINI* .
C**
C        ARGUMENT : LDFRAN (Entree) ==> Vrai/Faux pour francais/anglais.
C                           Input       True/False for french/english
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER INIMES, IREP, INUMER
C
      LOGICAL LDFRAN
#include "lficom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIFRA_MT',0,ZHOOK_HANDLE)
      IF (LFI%LFIFRA_LLPREA) THEN
C
C        Au premier appel du sous-programme, on initialise les variables
C     globales du logiciel (si cela n'a pas deja ete fait) .
C
C        At first routine call, initialisation of software global
C     variables (if not already done) .
C
        CALL LFIINI_MT (LFI, 2)
        LFI%LFIFRA_LLPREA=.FALSE.
      ENDIF
C
      LFI%LFRANC=LDFRAN
C
C        MESSAGERIE EVENTUELLE . MESSAGE, IF NECESSARY .
C
      IF (LFI%NIMESG.NE.2)  THEN 
        IF (LHOOK) CALL DR_HOOK('LFIFRA_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      INUMER=LFI%JPNIL
      INIMES=2
      IREP=0
      CLNSPR='LFIFRA'
      WRITE (UNIT=CLMESS,FMT='(''LDFRAN= '',L1)') LDFRAN
      CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,.FALSE.,CLMESS,
     S                CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFIFRA_MT',1,ZHOOK_HANDLE)
      END

