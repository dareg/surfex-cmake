C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE FATRAN(KREP,  KNUMER,  PCHAME, PCHAMS, LDOPT )
      USE FA_MOD, ONLY : FA_COM_DEFAULT,
     S                   FA_COM_DEFAULT_INIT,
     S                   NEW_FA_DEFAULT
      IMPLICIT NONE
#include "precision.h"
      INTEGER KREP, KNUMER
      LOGICAL LDOPT
      REAL (KIND=JPDBLR) PCHAME(*), PCHAMS(*)
      IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()
      CALL FATRAN_MT (FA_COM_DEFAULT, KREP,  KNUMER,  PCHAME, 
     S                PCHAMS, LDOPT )
      END SUBROUTINE

