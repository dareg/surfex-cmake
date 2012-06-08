C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE FAIENC(KREP, KNUMER, CDPREF, KNIVAU, CDSUFF,
     S                  PCHAMP, LDCOSP )
      USE FA_MOD, ONLY : FA_COM_DEFAULT,
     S                   FA_COM_DEFAULT_INIT,
     S                   NEW_FA_DEFAULT
      IMPLICIT NONE
#include "precision.h"
      INTEGER KREP, KNUMER, KNIVAU
      REAL (KIND=JPDBLR) PCHAMP (*)
      CHARACTER CDPREF*(*), CDSUFF*(*)
      LOGICAL LDCOSP
      IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()
      CALL FAIENC_MT (FA_COM_DEFAULT, KREP, KNUMER, CDPREF, KNIVAU, 
     S                CDSUFF, PCHAMP, LDCOSP )
      END SUBROUTINE

