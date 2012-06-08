C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE FACDEC(KREP, PA, PMIN, KNBIT, KDEC)
      USE FA_MOD, ONLY : FA_COM_DEFAULT,
     S                   FA_COM_DEFAULT_INIT,
     S                   NEW_FA_DEFAULT
      IMPLICIT NONE
#include "precision.h"
      REAL (KIND=JPDBLR) PA, PMIN
      INTEGER KREP, KNBIT, KDEC
      IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()
      CALL FACDEC_MT (FA_COM_DEFAULT, KREP, PA, PMIN, KNBIT, KDEC)
      END SUBROUTINE

