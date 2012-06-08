C Jan-2011 P. Marguinaud Interface to thread-safe LFI
      SUBROUTINE LFIDAH ( KDATE, KHEURE )
      USE LFIMOD, ONLY : LFICOM_DEFAULT,
     S                   LFICOM_DEFAULT_INIT,
     S                   NEW_LFI_DEFAULT
      IMPLICIT NONE
#ifndef f77
#include "precision.h"
      INTEGER (KIND=JPDBLE) KDATE, KHEURE
#else
      INTEGER KDATE, KHEURE
#endif
      IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()
      CALL LFIDAH_MT (LFICOM_DEFAULT, KDATE, KHEURE )
      END SUBROUTINE

