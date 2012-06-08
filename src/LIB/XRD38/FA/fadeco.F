C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE FADECO(KREP,   KNUMER, CDPREF, KNIVAU, CDSUFF,
     S                  LDCOSP, CDNOMA, KLNOMA, KVALCO, KLONGD,
     S                  PCHAMP )
      USE FA_MOD, ONLY : FA_COM_DEFAULT,
     S                   FA_COM_DEFAULT_INIT,
     S                   NEW_FA_DEFAULT
      IMPLICIT NONE
#include "precision.h"
      INTEGER KREP, KNUMER, KNIVAU, KLNOMA, KLONGD
      REAL (KIND=JPDBLR) PCHAMP (*)
      INTEGER (KIND=JPDBLE) KVALCO(*)
      LOGICAL LDCOSP
      CHARACTER CDPREF*(*), CDSUFF*(*), CDNOMA*(*)
      IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()
      CALL FADECO_MT (FA_COM_DEFAULT, KREP,   KNUMER, CDPREF, KNIVAU, 
     S                CDSUFF, LDCOSP, CDNOMA, KLNOMA, KVALCO, KLONGD,
     S                PCHAMP )
      END SUBROUTINE

