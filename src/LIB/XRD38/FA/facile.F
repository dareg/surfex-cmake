C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE FACILE(KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, 
     S                  KCHAMP, LDCOSP )
      USE FA_MOD, ONLY : FA_COM_DEFAULT,
     S                   FA_COM_DEFAULT_INIT,
     S                   NEW_FA_DEFAULT
      IMPLICIT NONE
#include "precision.h"
      INTEGER KREP, KNUMER, KNIVAU
      INTEGER (KIND=JPDBLE) KCHAMP (*)
      LOGICAL LDCOSP
      CHARACTER CDPREF*(*), CDSUFF*(*)
      IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()
      CALL FACILE_MT (FA_COM_DEFAULT, KREP, KNUMER, CDPREF, KNIVAU, 
     S                CDSUFF, KCHAMP, LDCOSP )
      END SUBROUTINE

