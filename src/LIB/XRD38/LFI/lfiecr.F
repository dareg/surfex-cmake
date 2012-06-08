C Jan-2011 P. Marguinaud Interface to thread-safe LFI
      SUBROUTINE LFIECR ( KREP, KNUMER, CDNOMA, KTAB, KLONG )
      USE LFIMOD, ONLY : LFICOM_DEFAULT,
     S                   LFICOM_DEFAULT_INIT,
     S                   NEW_LFI_DEFAULT
      IMPLICIT NONE
      CHARACTER CDNOMA*(*)
      INTEGER KREP, KNUMER, KLONG
#ifndef f77
#include "precision.h"
      INTEGER (KIND=JPDBLE)  KTAB (KLONG)
#else
      INTEGER KTAB (KLONG)
#endif
      IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()
      CALL LFIECR_MT (LFICOM_DEFAULT, KREP, KNUMER, 
     S                CDNOMA, KTAB, KLONG )
      END SUBROUTINE

