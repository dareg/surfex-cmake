C Jan-2011 P. Marguinaud Interface to thread-safe LFI
      SUBROUTINE LFILED ( KREP, KRANG, KTAB, KLONG, KRGPIM, 
     S                      KPOSEX, KRETIN )
      USE LFIMOD, ONLY : LFICOM_DEFAULT,
     S                   LFICOM_DEFAULT_INIT,
     S                   NEW_LFI_DEFAULT
      IMPLICIT NONE
      INTEGER KREP, KRANG, KLONG, KRGPIM, KPOSEX, KRETIN
#ifndef f77
#include "precision.h"
      INTEGER (KIND=JPDBLE)  KTAB (KLONG)
#else
      INTEGER KTAB (KLONG)
#endif
      IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()
      CALL LFILED_MT (LFICOM_DEFAULT, KREP, KRANG, KTAB, KLONG, KRGPIM, 
     S                      KPOSEX, KRETIN )
      END SUBROUTINE

