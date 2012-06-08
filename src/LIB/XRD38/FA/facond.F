C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE FACOND(KREP,   KNUMER, CDPREF, KNIVAU, CDSUFF,
     S                  PCHAMP, LDCOSP, CDNOMA, KLNOMA, PVALCO,
     S                  KLONGD )
      USE FA_MOD, ONLY : FA_COM_DEFAULT,
     S                   FA_COM_DEFAULT_INIT,
     S                   NEW_FA_DEFAULT
      IMPLICIT NONE
#include "precision.h"
      INTEGER KREP, KNUMER, KNIVAU, KLNOMA, KLONGD
      REAL (KIND=JPDBLR) PCHAMP (*), PVALCO (*)
      CHARACTER CDPREF*(*), CDSUFF*(*), CDNOMA*(*)
      LOGICAL LDCOSP
      IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()
      CALL FACOND_MT (FA_COM_DEFAULT, KREP,   KNUMER, CDPREF, KNIVAU, 
     S                CDSUFF, PCHAMP, LDCOSP, CDNOMA, KLNOMA, PVALCO,
     S                KLONGD )
      END SUBROUTINE

