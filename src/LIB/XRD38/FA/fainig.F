C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE FAINIG(KREP, KRANG,  CDPREF, KNIVAU, CDSUFF,
     S                  LDCOSP, KLCHAM, KSEC1, KSEC2, PSEC2, KSEC3,
     S                  PSEC3, KSEC4 )
      USE FA_MOD, ONLY : FA_COM_DEFAULT,
     S                   FA_COM_DEFAULT_INIT,
     S                   NEW_FA_DEFAULT
      IMPLICIT NONE
#include "precision.h"
      REAL (KIND=JPDBLR) PSEC3(*), PSEC2(*)
      INTEGER KREP, KRANG, KNIVAU, KLCHAM, KSEC1(*), KSEC2(*), KSEC3(*)
      INTEGER KSEC4(*)
      CHARACTER CDPREF*(*), CDSUFF*(*)
      LOGICAL LDCOSP
      IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()
      CALL FAINIG_MT (FA_COM_DEFAULT, KREP,   KRANG,  CDPREF, KNIVAU, 
     S                CDSUFF, LDCOSP, KLCHAM, KSEC1, KSEC2, PSEC2, 
     S                KSEC3, PSEC3, KSEC4 )
      END SUBROUTINE

