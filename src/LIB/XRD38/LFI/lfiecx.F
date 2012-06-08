C Jan-2011 P. Marguinaud Interface to thread-safe LFI
      SUBROUTINE LFIECX ( KREP, KRANG, KREC, KZONE, 
     S                      LDADON, KRETIN )
      USE LFIMOD, ONLY : LFICOM_DEFAULT,
     S                   LFICOM_DEFAULT_INIT,
     S                   NEW_LFI_DEFAULT
      IMPLICIT NONE
#ifndef f77
#include "precision.h"
      INTEGER (KIND=JPDBLE) KZONE (LFICOM_DEFAULT%JPLARX)
#else
      INTEGER KZONE (LFICOM_DEFAULT%JPLARX)
#endif
      INTEGER KREP, KRANG, KREC, KRETIN
      LOGICAL LDADON
      IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()
      CALL LFIECX_MT (LFICOM_DEFAULT, KREP, KRANG, KREC, KZONE, 
     S                      LDADON, KRETIN )
      END SUBROUTINE

