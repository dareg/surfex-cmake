C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIOFD_MT (LFI, KFACMD )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'obtenir la valeur courante
C     du FACTEUR MULTIPLICATIF par DEFAUT du logiciel LFI.
C**
C        ARGUMENT : KFACMD (Sortie) ==> Facteur multiplicatif par defaut
C
C        Le facteur multiplicatif est, pour une unite logique LFI, le
C     rapport entre la longueur d'enregistrement du fichier et LFI%JPRECL.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KFACMD, INIMES, IREP, INUMER
#include "lficom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIOFD_MT',0,ZHOOK_HANDLE)
      IF (LFI%LFIOFD_LLPREA) THEN
        CALL LFIINI_MT (LFI, 2)
        LFI%LFIOFD_LLPREA=.FALSE.
      ENDIF
C
      IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERGLA,'ON')
      KFACMD=LFI%MFACTU(0)
      IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERGLA,'OFF')
      INIMES=LFI%NIMESG
C
      IF (INIMES.EQ.2) THEN
        IREP=0
        INUMER=LFI%JPNIL
        CLNSPR='LFIOFD'
        WRITE (UNIT=CLMESS,FMT='(''KFACMD='',I3)') KFACMD
        CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,.FALSE.,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIOFD_MT',1,ZHOOK_HANDLE)
      END

