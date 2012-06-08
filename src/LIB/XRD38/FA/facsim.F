C Jan-2011 P. Marguinaud Interface to thread-safe FA
      SUBROUTINE FACSIM(KREP, KRANG, PCHAME, PCHAMS, 
     S                   KPULAS, KSTRON )
      USE FA_MOD, ONLY : FA_COM_DEFAULT,
     S                   FA_COM_DEFAULT_INIT,
     S                   NEW_FA_DEFAULT
      IMPLICIT NONE
#include "precision.h"
      INTEGER KREP, KRANG, KPULAS, KSTRON
      REAL (KIND=JPDBLR) PCHAME (*), PCHAMS (*)
      IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()
      CALL FACSIM_MT (FA_COM_DEFAULT, KREP, KRANG, PCHAME, PCHAMS, 
     S                KPULAS, KSTRON )
      END SUBROUTINE

