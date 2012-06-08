C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE FADECI(KREP,   KRANG,  CDNOMA, KVALCO, KLONGA,
     S                  KCHAMP, LDCOSP )
      USE FA_MOD, ONLY : FA_COM_DEFAULT,
     S                   FA_COM_DEFAULT_INIT,
     S                   NEW_FA_DEFAULT
      IMPLICIT NONE
#include "precision.h"
      INTEGER KREP, KRANG, KLONGA
      INTEGER (KIND=JPDBLE) KVALCO(*), KCHAMP(*)
      LOGICAL LDCOSP
      CHARACTER CDNOMA*(*)
      IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()
      CALL FADECI_MT (FA_COM_DEFAULT, KREP,   KRANG,  CDNOMA, KVALCO, 
     S                KLONGA, KCHAMP, LDCOSP )
      END SUBROUTINE

